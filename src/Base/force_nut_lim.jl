## ------------------------------------------------------------------
export _lb_biom_dual_prices
function _lb_biom_dual_prices(lep, 
        exch_idxs::Vector{Int}; 
        npoints = 5, 
        solver = Gurobi.Optimizer, 
        delta = 1e-2
    )

    lim_facs = Float64[]
    opm = FBAOpModel(lep, solver)
    for idx in exch_idxs
        v0 = lb(lep, idx)
        v1 = v0 * (1.0 - delta)
        # @show lep.colids[idx]
        @assert v0 < 0 && v1 < 0
        test_points = range(v0, v1; length = npoints)
        
        # R2FBAOpModel
        obj_m, obj_err, vars_ms, vars_errs = lb_dual_prices(opm, idx, test_points)
        push!(lim_facs, obj_m)
        # println()
        # println("-"^30)
        # @show colids(lep, idx)
        # @show test_points
        # @show obj_m
        # @show obj_err
    end
    return lim_facs
end

## ------------------------------------------------------------------
# For a fix biomass, which fluxes react to a fluctuation of a given target
export _sensible_fluxs
function _sensible_fluxs(net, flx_idx::Int, biom_idx::Int, biom0, biom1; 
        atol = 1e-3, npoints = 5, solver = Gurobi.Optimizer
    )
    opm = R2FBAOpModel(net, solver)
    # fix biom
    bounds!(opm, biom_idx, biom0, biom1)
    v0 = lb(net, flx_idx)
    test_points = range(v0, v0 * 0.9; length = npoints)
    obj_m, obj_err, vars_ms, vars_errs = flux_dual_prices(opm, flx_idx, test_points)
    return findall(abs.(vars_ms) .> atol)
end

## ------------------------------------------------------------------
# entropy driven contextualization
export _force_nut_limited
function _force_nut_limited(lep, glc_id, biom_id, exchs_ids; 
        biom_safe_factor = 0.8,
        ko_factor = 0.5, 
        protect_idxs = [], 
        oniter = (simdat) -> nothing,
        onsuccess = (simdat) -> nothing,
        niters = 800, 
        solver = Gurobi.Optimizer, 
        simdat = Dict()
    )

    lep = simdat["lep"] = deepcopy(lep)
    simdat["lep0"] = deepcopy(lep)
    simdat["glc_id"] = glc_id
    simdat["biom_id"] = biom_id
    simdat["exchs_ids"] = exchs_ids
    simdat["biom_safe_factor"] = biom_safe_factor
    simdat["ko_factor"] = ko_factor
    
    exchs_idxs = colindex.([lep], exchs_ids)
    simdat["exchs_idxs"] = exchs_idxs

    glc_idx = colindex(lep, glc_id)
    biom_idx = colindex(lep, biom_id)

    # max biom
    opm = FBAOpModel(lep, solver)
    optimize!(opm)
    biom1 = solution(opm, biom_idx)
    biom0 = biom1 * biom_safe_factor
    @show biom1

    simdat["biom0"] = biom0
    simdat["biom1"] = biom1

    # Protect biom
    bounds!(lep, biom_idx, biom0, biom1)

    # ---------------------------------
    # double prices
    sen_idxs = nothing
    lim_facs = nothing
    rand_idx = nothing
    m_glcs, max_ms = [], []
    traj_idxs, traj_b0s = [], []
    l0, u0 = nothing, nothing

    
    for it in 1:niters
        
        simdat["it"] = it
        
        @time try
            
            println("\n", "="^50)
            @show it
            
            oniter(simdat)
            
            lim_facs = _lb_biom_dual_prices(lep, exchs_idxs; solver, npoints = 5)
            simdat["lim_facs"] = lim_facs

            println("\n", "."^30)
            m_glc, max_m = 0.0, 0.0
            for (m, id) in zip(lim_facs, exchs_ids)
                if (id == glc_id); m_glc = abs(m);
                    elseif (abs(m) > 1e-3); nothing;
                    else; continue; 
                end
                max_m = max(max_m, abs(m))
                println(id, " :\t", round(m; digits = 3))
            end
            push!(m_glcs, m_glc)
            simdat["m_glcs"] = m_glcs
            push!(max_ms, max_m)
            simdat["max_ms"] = max_ms

            # success
            if (m_glc > 1e-3) && (max_m == m_glc)
                onsuccess(simdat)
                simdat["status"] = :success
                println("\n", "-"^30)
                @info("Limited :)")
                println()
                return simdat
            end

            println()
            sen_idxs = _sensible_fluxs(lep, glc_idx, biom_idx, biom0, biom1; solver, atol = 1e-3)
            sen_idxs = setdiff(sen_idxs, exchs_idxs, [biom_idx], protect_idxs)

            println("\n", "."^30)
            rand_idx = rand(sen_idxs)
            @show length(sen_idxs)
            @show lep.colids[rand_idx]

            l0, u0 = bounds(lep, rand_idx)
            @show ko_factor
            @show l0, u0
            bounds!(lep, rand_idx, l0 * ko_factor, u0 * ko_factor)
            l1, u1 = bounds(lep, rand_idx)
            @show l1, u1
            println()

            # Check biom
            opm = FBAOpModel(lep, solver)
            optimize!(opm)
            biom1 = solution(opm, biom_idx)

            println("\n", "."^30)
            @show biom0
            @show biom1
            simdat["biom0"] = biom0
            simdat["biom1"] = biom1

            # save simdat
            push!(traj_idxs, rand_idx)
            push!(traj_b0s, (l0, u0))
            simdat["traj_idxs"] = traj_idxs
            simdat["traj_b0s"] = traj_b0s

            println()
            
        catch err

            (err isa InterruptException) && rethrow(err)
            
            # ignore
            println("\n", "!"^30)
            @error err
            println()
            
            simdat["status"] = :error
            return simdat

        end # time

    end

    simdat["status"] = :unsuccess
    return simdat

end 
