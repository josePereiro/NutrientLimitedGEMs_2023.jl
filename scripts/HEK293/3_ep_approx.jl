## ------------------------------------------------------------------
@time begin
    using NutrientLimitedGEMs
    const NL = NutrientLimitedGEMs

    using ProjFlows
    using MetX
    using Gurobi
    using DataFrames
end

# ------------------------------------------------------------------
include("0_params.jl")
include("1_utils.jl")

# ------------------------------------------------------------------
let
    traj_dir = procdir(PROJ, ["HEK293", "sims"])
    frec = length(ARGS) > 0 ? parse(Int, ARGS[1]) : 5
    alg_ver = SIM_ID
    solver = LP_SOLVER
    recompute = false
    _compute_ep_data(traj_dir; solver, frec, alg_ver, recompute)
end

