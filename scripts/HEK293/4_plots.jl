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
# General 
let 
    traj_dir = procdir(NL, ["HEK293", "trajs"])
    p = plot()
    traj_lens, boxs = [], []
    biom1s, m_glcs = [], []
    @time for fn in readdir(traj_dir; join = true)
        traj = ldat(fn)
        traj["status"] == :success || continue
        get(traj, "dat_entropy_alg", "") == EP_ENTROPY_ALG_VERSION || continue

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
        title = "HEK293",
        xlabel = "ko steps", 
        ylabel = "max biom", 
        ylim = [0.0, 0.4]
    )
    push!(ps, p)

    p = scatter(traj_lens, m_glcs; 
        label = "", 
        title = "HEK293",
        xlabel = "ko steps", 
        ylabel = "glc m", 
        ylim = [0.0, 0.25],
    )
    push!(ps, p)

    p = histogram(m_glcs;
        bins = 80,
        label = "", 
        title = "HEK293",
        xlabel = "glc m", 
        ylabel = "count", 
    )
    push!(ps, p)

    ProjAssistant.sfig(NL, ps, 
        "HEK293", "trajectories", ".png"
    )

end

## ------------------------------------------------------------------
_colormap(xmin, xmax, x;)
## ------------------------------------------------------------------
## ------------------------------------------------------------------
# General 
let 
    ps = Plots.Plot[]
    traj_dir = procdir(NL, ["HEK293", "trajs"])
    p1 = plot(; xlabel = "ko steps", ylabel = "ΔS")
    Ssi1s, m_glc1s, ΔS1s = [], [], []
    ep_status1s = []
    @time for fn in readdir(traj_dir; join = true)
        traj = ldat(fn)
        traj["status"] == :success || continue
        get(traj, "dat_entropy_alg", "") == EP_ENTROPY_ALG_VERSION || continue

        # traj_idxs
        traj_idxs = traj["traj_idxs"]

        # Ss
        Ss = get(traj, "Ss", nothing)
        isnothing(Ss) && continue

        # Ss
        ep_status1 = last(traj["ep_statuses"])

        Ssis = findall(Ss) do S
            iszero(S) && return false
            isnan(S) && return false
            return true
        end
        isempty(Ssis) && continue
        @show Ssis

        # m_glc
        _m_glc = last(traj["m_glcs"])
        
        # ΔS
        ΔS = Ss[Ssis] .- first(Ss[Ssis])
        any(ΔS .> 15.0) && continue
        
        # Plots
        c = _colormap(0.0, 0.25, _m_glc; cname = "Grays")
        plot!(p1, Ssis, ΔS; 
            label = "", c, lw = 2, alpha = 0.2
        )

        push!(ΔS1s, last(ΔS))
        push!(Ssi1s, last(Ssis))
        push!(m_glc1s, _m_glc)
        push!(ep_status1s, ep_status1)

    end

    c = _colormap.(0.0, 0.25, m_glc1s; cname = "Grays")
    scatter!(p1, Ssi1s, ΔS1s; 
        label = "", c, m = 6, alpha = 0.8
    )
    push!(ps, p1)

    p2 = scatter(m_glc1s, ΔS1s; 
        label = "", c = :black, m = 6, alpha = 0.8,
        xlabel = "m glc", ylabel = "ΔS"
    )
    push!(ps, p2)

    ProjAssistant.sfig(NL, ps, 
        "HEK293", "entropy", ".png"
    )
end

## ------------------------------------------------------------------
## ------------------------------------------------------------------
let 
    traj_dir = procdir(NL, ["HEK293", "trajs"]) 
    p = plot()
    @time for fn in readdir(traj_dir; join = true)
        traj = ldat(fn)
        traj["status"] == :success || continue
        get(traj, "dat_entropy_alg", "") == EP_ENTROPY_ALG_VERSION || continue

        # Ss
        Ss = get(traj, "Ss", nothing)

        Ssis = findall(Ss) do S
            iszero(S) && return false
            isnan(S) && return false
            return true
        end
        Ss = Ss[Ssis]

        # ΔS
        ΔS = Ss .- first(Ss)
        any(ΔS .> 25.0) || continue

        # ep_status
        all_ok = all(traj["ep_statuses"] .== :converged)


        plot!(p, Ssis, ΔS; 
            label = "", 
            lw = 3, 
            color = all_ok ? :blue : :red
        )

    end 

    p
end