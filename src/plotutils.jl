# ------------------------------------------------------------------
export _colormap
function _colormap(p::Vector{RGB{Float64}}, xmin, xmax, x)
    N = length(p)
    r = range(xmin, xmax; length = N)
    _, idx = findmin((y) -> abs(x - y), r)
    return p[idx]
end

function _colormap(xmin, xmax, x;
        cname = "Grays", 
        N = 1000, 
        mid=0.5, 
        logscale=false,
    )
    cm = Plots.colormap(cname, N; mid, logscale)
    i = max(1, div(N, 5)):N
    _colormap(cm[i], xmin, xmax, x)
end

## ------------------------------------------------------------------
# Trajectoreis
export _plot_bash1 
function _plot_bash1(netid, ep_alg_version;
        biom_lims = (0.0, 0.7), 
        m_glcs_lims = (0.0, 0.10)
    )
    traj_dir = procdir(NutrientLimitedGEMs, [netid, "trajs"])
    p = plot()
    traj_lens, boxs = [], []
    biom1s, m_glcs = [], []
    @time for fn in readdir(traj_dir; join = true)
        endswith(fn, ".jls") || continue
        traj = ldat(fn)
        traj["status"] == :success || continue
        haskey(traj, ep_alg_version) || continue

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

    c = _colormap.(biom_lims..., biom1s; cname = "Grays")
    scatter_args = (;label = "", msc = :auto, title = netid, ms = 8, c)

    boxs ./= maximum(boxs)
    p = scatter(traj_lens, log10.(boxs); 
        scatter_args...,
        xlabel = "ko steps", 
        ylabel = "propto box vol",
        xlim = (0, Inf),
    )
    push!(ps, p)
    
    p = scatter(traj_lens, biom1s; 
        scatter_args...,
        xlabel = "ko steps", 
        ylabel = "max biom", 
        xlim = (0, Inf),
        ylim = biom_lims
    )
    push!(ps, p)

    p = scatter(traj_lens, m_glcs; 
        scatter_args...,
        xlabel = "ko steps", 
        ylabel = "glc m", 
        xlim = (0, Inf),
        ylim = m_glcs_lims,
    )
    push!(ps, p)

    p = histogram(m_glcs;
        bins = 80,
        label = "", c = :black,
        title = netid,
        xlabel = "glc m", 
        ylabel = "count", 
    )
    push!(ps, p)

    p = scatter(biom1s, m_glcs; 
        scatter_args...,
        xlabel = "max biom", 
        ylabel = "glc m", 
        xlim = biom_lims,
        ylim = m_glcs_lims
    )

    push!(ps, p)
    p = scatter(biom1s, log10.(boxs); 
        scatter_args...,
        xlabel = "max biom", 
        ylabel = "propto box vol", 
        xlim = biom_lims
    )
    push!(ps, p)

    sfig(NutrientLimitedGEMs, ps, 
        netid, "trajectories", ".png"
    )

end

## ------------------------------------------------------------------
# Entropy
export _plot_bash2
function _plot_bash2(netid, ep_alg_version;
        biom_lims
    )

    ps = Plots.Plot[]
    traj_dir = procdir(NutrientLimitedGEMs, [netid, "trajs"])
    p1 = plot(; xlabel = "ko steps", ylabel = "ΔS")
    p2 = plot(; xlabel = "ko steps", ylabel = "log vol box")
    Ssi1s, biom1s, ΔS1s, ΔV1s = [], [], [], []
    ep_status1s = []
    @time for fn in readdir(traj_dir; join = true)
        endswith(fn, ".jls") || continue
        traj = ldat(fn)
        traj["status"] == :success || continue
        epdat = get(traj, ep_alg_version, nothing)
        isnothing(epdat) && continue

        # Ss & V
        Ss = epdat["Ss"]
        Vs = epdat["box_vols"]
        
        # ep_status1
        ep_statuses = epdat["ep_statuses"]
        ep_status1 = last(ep_statuses)
        # ep_status1 == :converged || continue
        all(ep_statuses .== :converged) || continue

        # traj_idxs
        traj_idxs = traj["traj_idxs"]
        isempty(traj_idxs) && continue

        Ssis = findall(Ss) do S
            iszero(S) && return false
            isnan(S) && return false
            return true
        end
        isempty(Ssis) && continue
        # @show Ssis

        # bioms
        biom1 = traj["biom1"]
        
        # ΔS
        ΔSs = Ss[Ssis] .- first(Ss[Ssis])
        ΔVs = Vs[Ssis] ./ maximum(Vs[Ssis])
        # any(ΔS .> 5.0) && continue
        
        # Plots
        c = _colormap(biom_lims..., biom1; cname = "Grays")
        plot!(p1, Ssis, ΔSs; 
            label = "", c, lw = 2, alpha = 0.4, 
            ls = all(ep_statuses .== :converged) ? :solid : :dot
        )
        plot!(p2, Ssis, log10.(ΔVs); 
            label = "", c, lw = 2, alpha = 0.4, 
            ls = all(ep_statuses .== :converged) ? :solid : :dot
        )

        push!(ΔS1s, last(ΔSs))
        push!(ΔV1s, last(ΔVs))
        push!(Ssi1s, Int(last(Ssis)))
        push!(biom1s, biom1)
        push!(ep_status1s, ep_status1)

    end

    c = _colormap.(biom_lims..., biom1s; cname = "Grays")
    scatter!(p1, Ssi1s, ΔS1s; 
        label = "", c, m = 6, alpha = 0.8, 
        msc=:auto, xlim = (0, Inf)
    )
    push!(ps, p1)
    scatter!(p2, Ssi1s, log10.(ΔV1s); 
        label = "", c, m = 6, alpha = 0.8,
        msc=:auto, xlim = (0, Inf),
    )
    push!(ps, p2)

    p3 = scatter(biom1s, ΔS1s; 
        label = "", c, m = 6, alpha = 0.8,
        xlabel = "max biom", ylabel = "ΔS", 
        msc=:auto, xlim = biom_lims
    )
    push!(ps, p3)
    
    p4 = scatter(biom1s, log.(ΔV1s); 
        label = "", c, m = 6, alpha = 0.8,
        xlabel = "max biom", ylabel = "log vol box", 
        msc=:auto, xlim = biom_lims
    )
    push!(ps, p4)

    p5 = scatter(ΔS1s, log.(ΔV1s); 
        label = "", c, m = 6, alpha = 0.8,
        xlabel = "ΔS", ylabel = "log vol box", 
        msc=:auto
    )
    push!(ps, p5)

    p6 = scatter(Ssi1s, biom1s; 
        label = "", c, m = 6, alpha = 0.8,
        xlabel = "ok steps", ylabel = "max biom", 
        msc=:auto, 
        xlim = (0, Inf),
        ylim = biom_lims
    )
    push!(ps, p6)

    # write
    sfig(NutrientLimitedGEMs, ps, netid, "entropy", ".png")
end
