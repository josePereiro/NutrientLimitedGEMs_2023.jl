## ------------------------------------------------------------------
@time begin
    using NutrientLimitedGEMs
    const NL = NutrientLimitedGEMs
    
    using ProjFlows
    using MetX
    using Gurobi
end

# ------------------------------------------------------------------
include("0_params.jl")
include("1_utils.jl")

## ------------------------------------------------------------------
# Compute trajectories
let
    global lep0 = _setup_ecoli()
    exchs_ids0 = ["EX_GLC", "EX_NH4"]
    exchs_ids = extras.([lep0], exchs_ids0)
    exchs_idxs = colindex.([lep0], exchs_ids)
    
    protect_ids = filter(lep0.colids) do id
        startswith(id, "BIOMASS_Ecoli_core_w_GAM") && return true
        startswith(id, "EX_") && return true
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
        sim_hash = hash(sim["traj_idxs"])
        fn = procdir(PROJ, 
            ["ECOLI_CORE", "sims"], 
            (;hash = sim_hash), 
            ".sim.jls"
        )
        !isfile(fn) && sdat(sim, fn; verbose = true)
        println()

    end
    
end

## ------------------------------------------------------------------