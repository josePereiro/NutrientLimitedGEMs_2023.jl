## ------------------------------------------------------------------
@time begin
    using NutrientLimitedGEMs
    const NL = NutrientLimitedGEMs

    using ProjAssistant
    using MetXBase
    
    using Plots
    
end

# ------------------------------------------------------------------
include("0_params.jl")

## ------------------------------------------------------------------
let 
    fn = "/Users/Pereiro/.julia/dev/NutrientLimitedGEMs/data/processed/ECOLI/trajs/7659004865489147521.jls"
    global traj = ldat(fn)

end
## ------------------------------------------------------------------
# General 
let 
    traj_dir = procdir(NL, ["ECOLI", "trajs"])
    p = plot()
    traj_lens, boxs = [], []
    biom1s, m_glcs = [], []
    @time for fn in readdir(traj_dir; join = true)
        endswith(fn, ".jls") || continue
        @show fn
        traj = ldat(fn)
        traj["status"] == :success || continue
        # get(traj, "dat_entropy_alg", "") == EP_ENTROPY_ALG_VERSION || continue

        # traj_lens
        traj_idxs = traj["traj_idxs"]
        push!(traj_lens, length(traj_idxs))
        
        # box vol
        net = traj["net"]
        vol = prod(big.(net.ub .- net.lb))
        push!(boxs, vol)
        
        # biom
        biom1 = traj["biom1"]
        push!(biom1s, biom1)
        
        # glc_m
        _m_glcs = traj["m_glcs"]
        push!(m_glcs, last(_m_glcs))

    end

    ps = Plots.Plot[]

    boxs ./= maximum(boxs)
    p = scatter(traj_lens, log10.(boxs); 
        label = "", 
        xlabel = "ko steps", 
        ylabel = "propto box vol"
    )
    push!(ps, p)
    
    p = scatter(traj_lens, biom1s; 
        label = "", 
        title = "ECOLI",
        xlabel = "ko steps", 
        ylabel = "max biom", 
        ylim = [0.0, 0.4]
    )
    push!(ps, p)

    p = scatter(traj_lens, m_glcs; 
        label = "", 
        title = "ECOLI",
        xlabel = "ko steps", 
        ylabel = "glc m", 
        ylim = [0.0, 0.25],
    )
    push!(ps, p)

    p = histogram(m_glcs;
        bins = 80,
        label = "", 
        title = "ECOLI",
        xlabel = "glc m", 
        ylabel = "count", 
    )
    push!(ps, p)

    ProjAssistant.sfig(NL, ps, 
        "ECOLI", "trajectories", ".png"
    )

end