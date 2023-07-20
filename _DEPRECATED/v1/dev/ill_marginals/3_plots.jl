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
# Find unvalid trajs
illtraj_fns = _collect_ill_trajs(TRAJ_DIR, EP_ENTROPY_ALG_VERSION, ILL_ΔSTH)

## ------------------------------------------------------------------
# Marginals
let 

    for fn in sort(illtraj_fns)

        println("\n", "="^60)
        @show fn
        traj = ldat(fn)
        net = traj["net"]
        @show size(net)

        alg = basename(fn)
        _plot_hists(alg)
    end

end