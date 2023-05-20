## ------------------------------------------------------------------
@time begin
    using NutrientLimitedGEMs
    const NL = NutrientLimitedGEMs
    
    using ProjFlows
    using ContextDBs
    using MetX
    using Gurobi
end

## ------------------------------------------------------------------
include("0_params.jl")
include("1_utils.jl")

## ------------------------------------------------------------------
# Compute trajectories
# @save let # TODO: make this happend. Capture all assigments

# TODO: Mantra: Just like comments, the input workflow should not interfer with an existing program
# TODO: names: save! -> track! -> commit
# IDEA: Use definition diff for tracking track defvars0 - defvars1 since a given point. eg. a context change
# Git commit every time the top context change/run
let
    @context "KO_MC"

    @save! lep0 = _setup_ecoli()
    @save! exchs_ids0 = ["EX_GLC", "EX_NH4"]
    @save! exchs_ids = extras.([lep0], exchs_ids0)
    @save! exchs_idxs = colindex.([lep0], exchs_ids)
    
    @save! protect_ids = filter(lep0.colids) do id
        startswith(id, "BIOMASS_Ecoli_core_w_GAM") && return true
        startswith(id, "EX_") && return true
        return false
    end
    @save! protect_idxs = colindex.([lep0], protect_ids)
    
    glc_id, biom_id = extras.([lep0], ["EX_GLC", "BIOM"])
    
    # ---------------------------------
    @save! glc_id biom_id
    @save! biom_safe_factor = 0.1
    @save! niters = 500
    @save! ko_factor = KO_FACTORs
    @save! solver = LP_SOLVER

    # TODO: finish this
    # Think how to interact with methods
    
    while true

        sim = _force_nut_limited(lep0, 
            glc_id, biom_id, exchs_ids; 
            biom_safe_factor,
            protect_idxs,
            ko_factor,
            niters,
            solver
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