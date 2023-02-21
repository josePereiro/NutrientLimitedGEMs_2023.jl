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

# ------------------------------------------------------------------
include("0_params.jl")
include("1_utils.jl")

## ------------------------------------------------------------------
let
    traj_dir = procdir(NL, ["HEK293", "trajs"])
    frec = length(ARGS) > 0 ? parse(Int, ARGS[1]) : 5
    alg_ver = EP_ENTROPY_ALG_VERSION
    solver = LP_SOLVER
    _compute_ep_entropy_data(traj_dir; solver, frec, alg_ver)
end