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
    
    net0 = _setup_heknet()
    
    exchs_ids0 = [
        "EX_GLC", "EX_NH4", "EX_TYR","EX_HIS","EX_ASN","EX_ARG",
        "EX_ILE","EX_PHE","EX_MET","EX_THR","EX_LYS","EX_LEU",
        "EX_VAL","EX_TRP","EX_GLU","EX_GLN"
    ]
    exchs_ids = extras.([net0], exchs_ids0)
    exchs_idxs = rxnindex.([net0], exchs_ids)

    protect_ids = filter(net0.rxns) do id
        startswith(id, "R_biomass") && return true
        startswith(id, "R_EX_") && return true
        return false
    end
    protect_idxs = rxnindex.([net0], protect_ids)

    glc_id = extras(net0, "EX_GLC")
    biom_id = extras(net0, "BIOM")

    # ---------------------------------
    while true

        traj = _force_nut_limited(net0, 
            glc_id, biom_id, exchs_ids; 
            biom_safe_factor = 0.1,
            protect_idxs,
            ko_factor = KO_FACTOR,
            niters = 500,
            solver = LP_SOLVER
        )
        
        traj["status"] == :success || continue
        traj_hash = hash(traj["traj_idxs"])
        fn = procdir(NL, ["HEK293", "trajs"], traj_hash, ".jls")
        !isfile(fn) && sdat(traj, fn)
    end
    
end

## ------------------------------------------------------------------


# ------------------------------------------------------------------
