## ------------------------------------------------------------------
@time begin
    using NutrientLimitedGEMs_2023
    const NL = NutrientLimitedGEMs_2023

    using ProjFlows
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
    traj_dir = procdir(NL, ["IJR904", "sims"])
    frec = length(ARGS) > 0 ? parse(Int, ARGS[1]) : 5
    alg_ver = SIM_ID
    solver = LP_SOLVER
    recompute = false
    epsconv = 1e-5
    _compute_ep_data(traj_dir; solver, frec, alg_ver, recompute, epsconv)
end