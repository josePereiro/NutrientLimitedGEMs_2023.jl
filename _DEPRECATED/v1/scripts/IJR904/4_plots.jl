## ------------------------------------------------------------------
@time begin
    using NutrientLimitedGEMs_2023

    using ProjFlows
    using Gurobi
    using MetXBase
    
    using Plots
    
end

## ------------------------------------------------------------------
include("0_params.jl")

# ------------------------------------------------------------------
# Trajectories
_plot_bash1("IJR904", SIM_ID;
    biom_lims = (0.0, 0.4), 
    m_glcs_lims = (0.0, 0.30)
)

## ------------------------------------------------------------------
# Entropy
_plot_bash2("IJR904", SIM_ID;
    biom_lims = (0.0, 0.4), 
)
