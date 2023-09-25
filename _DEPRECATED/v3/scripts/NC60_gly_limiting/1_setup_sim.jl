@time begin
    using NutrientLimitedGEMs
    using ProjFlows
    using Gurobi
    using ContextDBs
    using Base.Threads
end

## ------------------------------------------------------------------
# PROJECT
## ------------------------------------------------------------------

PROJ = Project0(NutrientLimitedGEMs)

## ------------------------------------------------------------------
# ContextDB
## ------------------------------------------------------------------

@newcontextdb!
cacherefs_dir!(cachedir(PROJ))

## ------------------------------------------------------------------
# TOP context

@context! WIP = v"0.1.0"

#= DESCRIPTION:

1. Create the ko bundles starting from HumanGEM1. 
Each bondle will check the biomass and do fva on the nen_essential intakes. 

=#

## ------------------------------------------------------------------
@tempcontext ["GLOBALS"] begin
    @stage! LP_SOLVER = Gurobi.Optimizer
    @stage! NTHREADS = clamp(nthreads() - 1, 1, 2)
end


return nothing