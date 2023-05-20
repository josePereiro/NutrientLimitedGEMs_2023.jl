## ------------------------------------------------------------------
# Compute EP stuff
export _compute_ep_data
function _compute_ep_data(traj_dir;
        solver = Gurobi.Optimizer, 
        frec = 5, 
        alg_ver = "EP", 
        niters = 15000, 
        epsconv = 1e-8,
        recompute = false, 
    )

    done = Set{String}()
    for _ in 1:niters
        
        # traj file
        files = readdir(traj_dir; join = true)
        isempty(files) && break
        length(files) == length(done) && break
        fn = rand(files)
        fn ∈ done && continue
        push!(done, fn)
        endswith(fn, ".jls") || continue

        _, sim = ldat(fn)
        sim["status"] == :success || continue
        
        println("\n", "="^30)

        ko_factor = sim["ko_factor"]
        @show ko_factor
        traj_idxs = sim["traj_idxs"]
        L = length(traj_idxs) + 1
        @show L
        lep0 = deepcopy(sim["lep0"])
        
        dowrite = false
        
        # EP data
        epdat = get!(sim, alg_ver, Dict()) 
        Ss = get!(epdat, "Ss", zeros(L)) 
        Fs = get!(epdat, "Fs", zeros(L)) 
        log_ZQs = get!(epdat, "log_ZQs", zeros(L)) 
        ∑logZ_Qns = get!(epdat, "∑logZ_Qns", zeros(L)) 
        rxns_counts = get!(epdat, "rxns_counts", zeros(Int, L))
        ep_statuses = get!(epdat, "ep_statuses", fill(:unset, L))
        box_vols = get!(epdat, "box_vols", zeros(BigFloat, L))
        ep_avs = get!(epdat, "ep_avs", Vector{Vector{Float64}}(undef, L))
        ep_vas = get!(epdat, "ep_vas", Vector{Vector{Float64}}(undef, L))
        box_lbs = get!(epdat, "box_lbs", Vector{Vector{Float64}}(undef, L))
        box_ubs = get!(epdat, "box_ubs", Vector{Vector{Float64}}(undef, L))

        # compute entropy
        _compute_entropy = (i) -> let 
            
            lep = box(lep0, solver; 
                verbose = true, 
                reduce = true,
                eps = 1e-5 # TODO: Test smaller
            )
            @show size(lep)
            epm = FluxEPModelT0(lep)
            config!(epm; verbose = true, epsconv)
            converge!(epm)
            ep_status = convergence_status(epm)
            @show ep_status
            S = entropy(epm)
            @show S
            F, log_ZQ, ∑logZ_Qn = free_energy(epm)
            @show F

            # Up state
            Ss[i] = S
            Fs[i] = F
            log_ZQs[i] = log_ZQ
            ∑logZ_Qns[i] = ∑logZ_Qn
            rxns_counts[i] = size(lep, 2)
            ep_statuses[i] = ep_status
            box_vols[i] = prod(big.(lep.ub .- lep.lb))
            ep_avs[i] = mean(epm)
            ep_vas[i] = var(epm)
            box_lbs[i] = lb(lep)
            box_ubs[i] = ub(lep)
        end

        # entropy lep0
        try
            i = 0
            println("\n", "."^30)
            @show i
            @show frec

            _compute_entropy(1)
            
        catch err
            (err isa InterruptException) && rethrow(err)
            println("\n", "!"^30)
            @error err
            println()
        end


        # entropy sim
        for (i, idx) in enumerate(traj_idxs)
            try
                println("\n", "."^30)
                @show i
                @show frec
                
                # add ko
                l0, u0 = bounds(lep0, idx)
                @show colids(lep0, idx)
                @show (l0, u0)
                bounds!(lep0, idx, l0 * ko_factor, u0 * ko_factor)
                l1, u1 = bounds(lep0, idx)
                @show (l1, u1)

                # skip
                skip = !recompute && (Ss[i + 1] != 0.0) # Already computed
                skip |= !(i == 1 || (i + 1) == L || iszero(rem(i, frec)))
                skip && continue

                _compute_entropy(i + 1)
                
            catch err
                (err isa InterruptException) && rethrow(err)
                println("\n", "!"^30)
                # @error err
                rethrow(err)
                println()
            end
            dowrite = true

        end # for (i, idx)
        
        # save
        if dowrite 
            println("\n", "-"^30)
            sdat(sim, fn; verbose = true)
            println()
        end

    end # for sim in trajs
end