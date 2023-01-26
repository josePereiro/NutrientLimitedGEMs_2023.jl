## ------------------------------------------------------------------
@time begin
    using NutrientLimitedGEMs
    const NL = NutrientLimitedGEMs

    using ProjAssistant
    using MetXBase
    using MetXOptim
    using MetXNetHub
    using MetXEP
    using MetXCultureHub
    using Gurobi
    using Serialization
    using ArgParse
    
    using Plots
    using MetXPlots
    
    # Pkg.add("https://github.com/josePereiro/ImgTools.jl")
    using ImgTools
end

## ------------------------------------------------------------------
include("0_params.jl")
include("1_utils.jl")

## ------------------------------------------------------------------
# Compute Entropy
let

    frec = length(ARGS) > 0 ? parse(Int, ARGS[1]) : 5
    
    traj_dir = procdir(NL, ["HEK293", "trajs"])

    for _ in 1:5000
    
        fn = rand(readdir(traj_dir; join = true))

        global traj = ldat(fn)
        traj["status"] == :success || continue
        
        println("\n", "="^30)

        ko_factor = traj["ko_factor"]
        @show ko_factor
        traj_idxs = traj["traj_idxs"]
        L = length(traj_idxs)
        @show L
        net0 = deepcopy(traj["net0"])
        
        dowrite = false
        
        # EP data
        Ss = get!(traj, "Ss", zeros(L)) 
        rxns_counts = get!(traj, "rxns_counts", zeros(Int, L))
        ep_statuses = get!(traj, "ep_statuses", fill(:unset, L))

        # Reset if new version
        dat_entropy_alg = get(traj, "dat_entropy_alg", "")

        if (dat_entropy_alg != EP_ENTROPY_ALG_VERSION) 
            Ss .= 0.0
            rxns_counts .= 0
            ep_statuses .= :unset
            dowrite = true
        end

        for (i, idx) in enumerate(traj_idxs)

            try

                println("\n", "."^30)
                @show i
                @show frec

                # Already computed
                skip = (Ss[i] != 0.0)
                skip &= !(i == 1 || i == L || iszero(rem(i, frec)))
                skip && continue
                
                # Prepare net
                l0, u0 = bounds(net0, idx)
                @show net0.rxns[i]
                @show (l0, u0)
                bounds!(net0, idx, l0 * ko_factor, u0 * ko_factor)
                l1, u1 = bounds(net0, idx)
                @show (l1, u1)

                net = box(net0, LP_SOLVER; 
                    verbose = true, 
                    reduce = true,
                    protect_obj = true,
                    eps = 1e-5
                )
                @show size(net)
                epm = FluxEPModelT0(net)
                converge!(epm; verbose = true)
                ep_status = convergence_status(epm)
                @show ep_status
                S = entropy(epm)
                @show S
                
                # Up state
                Ss[i] = S
                rxns_counts[i] = size(net, 2)
                ep_statuses[i] = ep_status
                
            catch err

                (err isa InterruptException) && rethrow(err)

                println("\n", "!"^30)
                @error err
                println()
                
            end

            dowrite = true
        
        end # for (i, idx)

        # alg label
        traj["dat_entropy_alg"] = EP_ENTROPY_ALG_VERSION
        
        # save
        dowrite && sdat(traj, fn)

    end # for traj in trajs

end