# TODO: reproduce the analysis done using ecoli_core (see DEPRECATED)
## ------------------------------------------------------------------
LP_SOLVER = Gurobi.Optimizer
KO_FACTOR = 1e-2

## ------------------------------------------------------------------
TRAJ_DIR = procdir(NutrientLimitedGEMs, ["ECOLI_CORE", "trajs"])
EP_ENTROPY_ALG_VERSION = "EPv6"
ILL_ΔSTH = 20

## ------------------------------------------------------------------
return nothing