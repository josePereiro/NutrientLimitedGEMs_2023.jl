## ------------------------------------------------------------------
@time begin
    using NutrientLimitedGEMs
    using ProjFlows
    using Gurobi
    using MetXBase
    using MetX
    using Plots
end

# ------------------------------------------------------------------
include("0_params.jl")

## ------------------------------------------------------------------
# Trajectories
_plot_bash1(PROJ, "HEK293", SIM_ID; 
    biom_lims = (0.0, 0.4), 
    m_glcs_lims = (0.0, 0.30), 
    biom1_th = 0.0
)

## ------------------------------------------------------------------
# Entropy
_plot_bash2(PROJ, "HEK293", SIM_ID;
    biom_lims = (0.0, 0.4), 
    biom1_th = 0.0,
)

