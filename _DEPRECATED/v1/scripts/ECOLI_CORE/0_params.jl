# TODO: reproduce the analysis done using ecoli_core (see DEPRECATED)

## ------------------------------------------------------------------
# Project
PROJ = Project0(NutrientLimitedGEMs)

## ------------------------------------------------------------------
# TODO: load ContextDB

## ------------------------------------------------------------------
@initcontext! SIM_ID = "SIM_28.3.2023.v1"

## ------------------------------------------------------------------
@withcontext! ["PARAMS"] begin
    @save! LP_SOLVER = Gurobi.Optimizer
    @save! KO_EPS = 1e-6
end


return nothing