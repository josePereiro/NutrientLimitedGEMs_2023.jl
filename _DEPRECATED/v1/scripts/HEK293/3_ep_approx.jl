## ------------------------------------------------------------------
@time begin
    using NutrientLimitedGEMs_2023
    const NL = NutrientLimitedGEMs_2023

    using ProjFlows
    using MetX
    using Gurobi
    using DataFrames
end

# ------------------------------------------------------------------
include("0_params.jl")
include("1_utils.jl")

## ------------------------------------------------------------------
let
    traj_dir = procdir(PROJ, ["HEK293", "sims"])
    frec = length(ARGS) > 0 ? parse(Int, ARGS[1]) : 5
    alg_ver = SIM_ID
    solver = LP_SOLVER
    recompute = true
    epsconv = 1e-6
    ko_eps = 1e-4
    _compute_ep_data(traj_dir; epsconv, solver, frec, alg_ver, recompute, ko_eps)
end

