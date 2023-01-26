## ------------------------------------------------------------------
export _lb_biom_dual_prices
function _lb_biom_dual_prices(net, 
        exch_idxs::Vector{Int}, biom_idx::Int; 
        npoints = 5, 
        solver = Gurobi.Optimizer
    )

    lim_facs = Float64[]
    opm = FBAFluxOpModel(net, solver)
    for idx in exch_idxs
        v0 = lb(net, idx)
        v1 = v0 * 0.9
        # @show net.rxns[idx]
        @assert v0 < 0 && v1 < 0
        test_points = range(v0, v1; length = npoints)
        # @show test_points
        
        # R2FBAFluxOpModel
        _ms, _errs = lb_dual_prices(opm, idx, test_points)
        m = getindex(_ms, biom_idx)
        push!(lim_facs, m)
    end
    return lim_facs
end

## ------------------------------------------------------------------
# For a fix biomass, which fluxes react to a fluctuation of a given target
export _sensible_fluxs
function _sensible_fluxs(net, flx_idx::Int, biom_idx::Int, biom0, biom1; 
        atol = 1e-3, npoints = 5, solver = Gurobi.Optimizer
    )
    opm = R2FBAFluxOpModel(net, solver)
    # fix biom
    bounds!(opm, biom_idx, biom0, biom1)
    v0 = lb(net, flx_idx)
    test_points = range(v0, v0 * 0.9; length = npoints)
    ms, errs = flux_dual_prices(opm, flx_idx, test_points)
    return findall(abs.(ms) .> atol)
end

## ------------------------------------------------------------------
# entropy driven contextualization
export _force_nut_limited
function _force_nut_limited(net, glc_id, biom_id, exchs_ids; 
        biom_safe_factor = 0.8,
        ko_factor = 0.5, 
        protect_idxs = [], 
        oniter = (state) -> nothing,
        onsuccess = (state) -> nothing,
        niters = 800, 
        solver = Gurobi.Optimizer
    )

    state = Dict()

    net = state["net"] = deepcopy(net)
    state["net0"] = deepcopy(net)
    state["glc_id"] = glc_id
    state["biom_id"] = biom_id
    state["exchs_ids"] = exchs_ids
    state["biom_safe_factor"] = biom_safe_factor
    state["ko_factor"] = ko_factor
    
    exchs_idxs = rxnindex.([net], exchs_ids)
    state["exchs_idxs"] = exchs_idxs

    glc_idx = rxnindex(net, glc_id)
    biom_idx = rxnindex(net, biom_id)

    # max biom
    opm = FBAFluxOpModel(net, solver)
    optimize!(opm)
    biom1 = solution(opm, biom_idx)
    biom0 = biom1 * biom_safe_factor
    @show biom1

    state["biom0"] = biom0
    state["biom1"] = biom1

    # Protect biom
    bounds!(net, biom_idx, biom0, biom1)

    # ---------------------------------
    # double prices
    sen_idxs = nothing
    lim_facs = nothing
    rand_idx = nothing
    m_glcs, max_ms = [], []
    traj_idxs, traj_b0s = [], []
    l0, u0 = nothing, nothing

    for it in 1:niters

        state["it"] = it

        @time try

            println("\n", "="^50)
            @show it

            oniter(state)

            lim_facs = _lb_biom_dual_prices(net, exchs_idxs, biom_idx; solver, npoints = 5)
            state["lim_facs"] = lim_facs
            
            println("\n", "."^30)
            m_glc, max_m = 0.0, 0.0
            for (m, id) in zip(lim_facs, exchs_ids)
                if id == glc_id; m_glc = abs(m)
                    elseif (abs(m) > 1e-3) 
                    else; continue; 
                end
                max_m = max(max_m, abs(m))
                println(id, " :\t", round(m; digits = 3))
            end
            push!(m_glcs, m_glc)
            state["m_glcs"] = m_glcs
            push!(max_ms, max_m)
            state["max_ms"] = max_ms
            # success
            if (m_glc > 1e-3) && (max_m == m_glc)
                onsuccess(state)
                state["status"] = :success
                @info("Limited!!!")
                return state
            end
            
            println()
            sen_idxs = _sensible_fluxs(net, glc_idx, biom_idx, biom1 * 0.99, biom1; solver, atol = 1e-3)
            sen_idxs = setdiff(sen_idxs, exchs_idxs, [biom_idx], protect_idxs)
            
            println("\n", "."^30)
            rand_idx = rand(sen_idxs)
            @show length(sen_idxs)
            @show net.rxns[rand_idx]

            l0, u0 = bounds(net, rand_idx)
            @show ko_factor
            @show l0, u0
            bounds!(net, rand_idx, l0 * ko_factor, u0 * ko_factor)
            l1, u1 = bounds(net, rand_idx)
            @show l1, u1
            println()

            # Check biom
            opm = FBAFluxOpModel(net, solver)
            optimize!(opm)
            biom1 = solution(opm, biom_idx)
            
            println("\n", "."^30)
            @show biom0
            @show biom1
            state["biom0"] = biom0
            state["biom1"] = biom1
            
            # save state
            push!(traj_idxs, rand_idx)
            push!(traj_b0s, (l0, u0))
            state["traj_idxs"] = traj_idxs
            state["traj_b0s"] = traj_b0s

            println()
            
        catch err
            
            (err isa InterruptException) && rethrow(err)
            
            println("\n", "!"^30)
            @error err
            println()

            state["status"] = :error
            return state

        end # time

        
    end

    state["status"] = :unsuccess
    return state
    
end 


# ------------------------------------------------------------------
function _colormap(p::Vector{RGB{Float64}}, xmin, xmax, x)
    N = length(p)
    r = range(xmin, xmax; length = N)
    _, idx = findmin((y) -> abs(x - y), r)
    return p[idx]
end

function _colormap(xmin, xmax, x;
        cname = "Grays", 
        N = 1000, 
        mid=0.5, 
        logscale=false,
    )
    cm = Plots.colormap(cname, N; mid, logscale)
    _colormap(cm, xmin, xmax, x)
end


## ------------------------------------------------------------------
nothing