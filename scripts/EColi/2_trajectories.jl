## ------------------------------------------------------------------
@time begin
    using NutrientLimitedGEMs
    const NL = NutrientLimitedGEMs
    
    using ProjAssistant
    using MetXBase, MetXOptim, MetXNetHub
    using MetXEP, MetXCultureHub
    using Gurobi
    
    using Plots
    using MetXPlots
    
    # Pkg.add("https://github.com/josePereiro/ImgTools.jl")
    using ImgTools
end

# ------------------------------------------------------------------
include("0_params.jl")
include("1_utils.jl")

## ------------------------------------------------------------------
# Compute trajectories
let
    global net0 = _setup_ecoli()
    
    exchs_ids0 = ["EX_GLC", "EX_NH4"]
    exchs_ids = extras.([net0], exchs_ids0)
    exchs_idxs = rxnindex.([net0], exchs_ids)

    protect_ids = filter(net0.rxns) do id
        startswith(id, "BIOMASS_Ecoli_core_w_GAM") && return true
        startswith(id, "EX_") && return true
        return false
    end
    protect_idxs = rxnindex.([net0], protect_ids)

    glc_id = extras(net0, "EX_GLC")
    biom_id = extras(net0, "BIOM")

    # ---------------------------------
    for it in 1:150000

        global traj = _force_nut_limited(net0, 
            glc_id, biom_id, exchs_ids; 
            biom_safe_factor = 0.1,
            protect_idxs,
            ko_factor = KO_FACTOR,
            niters = 500,
            solver = LP_SOLVER
        )
        
        # traj["status"] == :success || continue
        traj["status"] == :success && break
        traj_hash = hash(traj["status"])
        sdat(NL, traj,
            ["ECOLI", "trajs"], traj_hash, ".jls"
        )
    end
    
end

## ------------------------------------------------------------------