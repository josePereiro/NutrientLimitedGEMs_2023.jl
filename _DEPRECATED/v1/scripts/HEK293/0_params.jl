## ------------------------------------------------------------------
begin
    using NutrientLimitedGEMs
    using ProjFlows
    using ContextDBs
    using Gurobi
    using Random
end

## ------------------------------------------------------------------
Random.seed!(123)

## ------------------------------------------------------------------
SIM_ID = "SIM_11.4.2023.Calzadilla.v1"
LP_SOLVER = Gurobi.Optimizer
KO_EPS = 1e-7

## ------------------------------------------------------------------
# SIM_ID = "SIM_4.2023.TrueReduction"
# LP_SOLVER = Gurobi.Optimizer
# KO_EPS = 0.0

## ------------------------------------------------------------------
# Project
PROJ = Project0(NutrientLimitedGEMs)

nothing