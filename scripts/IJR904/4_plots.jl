## ------------------------------------------------------------------
@time begin
    using NutrientLimitedGEMs

    using ProjAssistant
    using Gurobi
    using MetXBase
    
    using Plots
    
end

## ------------------------------------------------------------------
include("0_params.jl")

# ------------------------------------------------------------------
# Trajectories
_plot_bash1("IJR904", EP_ENTROPY_ALG_VERSION;
    biom_lims = (0.0, 0.4), 
    m_glcs_lims = (0.0, 0.30)
)

# ------------------------------------------------------------------
# Entropy
_plot_bash2("IJR904", EP_ENTROPY_ALG_VERSION;
    biom_lims = (0.0, 0.4), 
)