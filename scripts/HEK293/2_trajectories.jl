## ------------------------------------------------------------------
@time begin
    using NutrientLimitedGEMs
    using ProjFlows
    using MetX
    using Gurobi
    using DataFrames
end

# ------------------------------------------------------------------
include("0_params.jl")
include("1_utils.jl")

# ------------------------------------------------------------------
# Compute trajectories
let
    
    lep0 = _base_heknet(PROJ)
    
    exchs_ids0 = [
        "EX_GLC", "EX_NH4", "EX_TYR","EX_HIS","EX_ASN","EX_ARG",
        "EX_ILE","EX_PHE","EX_MET","EX_THR","EX_LYS","EX_LEU",
        "EX_VAL","EX_TRP","EX_GLU","EX_GLN"
    ]
    exchs_ids = extras.([lep0], exchs_ids0)
    exchs_idxs = colindex.([lep0], exchs_ids)

    protect_ids = filter(colids(lep0)) do id
        startswith(id, "R_biomass") && return true
        startswith(id, "R_EX_") && return true
        return false
    end
    protect_idxs = colindex.([lep0], protect_ids)

    glc_id, biom_id = extras.([lep0], ["EX_GLC", "BIOM"])

    # ---------------------------------
    while true

        sim = _force_nut_limited(lep0, 
            glc_id, biom_id, exchs_ids; 
            biom_safe_factor = 0.1,
            protect_idxs,
            ko_factor = KO_FACTOR,
            niters = 500,
            solver = LP_SOLVER
        )
        
        sim["status"] == :success || continue
        fn = procdir(PROJ, 
            ["HEK293", "sims"], 
            (;hash = hash(sim["traj_idxs"])), 
            ".sim.jls"
        )
        !isfile(fn) && sdat(sim, fn; verbose = true)
        println()
    end
    
end

## ------------------------------------------------------------------


# ------------------------------------------------------------------
