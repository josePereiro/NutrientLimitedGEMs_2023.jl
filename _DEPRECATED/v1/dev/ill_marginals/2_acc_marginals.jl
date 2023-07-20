## ------------------------------------------------------------------
@time begin
    using NutrientLimitedGEMs
    using MetX
    using MetXBase
    using MetXBase: _histogram!
    using MetXOptim
    using MetXNetHub
    using MetXEP
    using MetXMC
    using MetXCultureHub
    using Gurobi
    using Serialization
    using Plots
    using ProjAssistant
end

# -------------------------------------------------------------------
include("0_params.jl")
include("1_tools.jl")

## ------------------------------------------------------------------
# Marginals
let 
    # Find unvalid trajs
    println("\n", "="^60)
    illtraj_fns = _collect_ill_trajs(TRAJ_DIR, EP_ENTROPY_ALG_VERSION, ILL_ΔSTH)

    while true
        for fn in sort(illtraj_fns)[1:15]

            println("\n", "="^60)
            @show fn
            traj = ldat(fn)
            net = traj["net"]
            @show size(net)

            alg = basename(fn)
            _accumulate_marginals(net, alg; 
                save_frec = 60.0, 
                tout = 5.0 * 60.0
            )
        end
    end

end