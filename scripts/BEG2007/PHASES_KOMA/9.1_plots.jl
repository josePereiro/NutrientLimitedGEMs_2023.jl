## --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
@time begin
    using Random
    using MetXGEMs
    using MetXBase
    using ProjFlows
    using Statistics
    using CairoMakie
    using BlobBatches
    using Base.Threads
    using Combinatorics
    using Base.Threads: Atomic
    using NutrientLimitedGEMs
end

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
include("1_setup.jl")
include("2_utils.jl")

## --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# biomass vs downset length histogram
let
    _simver = "ECOLI-CORE-BEG2007-PHASE_I-0.1.0"

    n0 = 0 # init file
    n1 = Inf # non-ignored file count
    cid = (@__FILE__, _simver, "biomass:core vs gem", n0, n1)
    lk = ReentrantLock()
    _, ret = withcachedat(PROJ, :get!, cid) do
        _gem_biomv, _core_biomv = Float64[], Float64[]
        _th_readdir(_simver; n0, n1, nthrs = 1) do bbi, bb
            haskey(bb["meta"], "core_biomass_fba.ver") || return :ignore
            haskey(bb["meta"], "gem_biomass_fba.ver") || return :ignore
            # load frame
            feasets_db = bb["core_feasets"]
            _th_gem_biomv, _th_core_biomv = Float64[], Float64[]
            for feasets_blob0 in feasets_db
                for (_fealen, feasets_blob1) in feasets_blob0
                    _core_biom = feasets_blob1["core_biomass_fba.biom"]
                    _gem_biom  = feasets_blob1["gem_biomass_fba.biom"]
                    isnan(_core_biom) && continue
                    isnan(_gem_biom) && continue
                    push!(_th_core_biomv, _core_biom)
                    push!(_th_gem_biomv, _gem_biom)
                end
            end # for feasets_blob0
            lock(lk) do
                push!(_core_biomv, _th_core_biomv...)
                push!(_gem_biomv, _th_gem_biomv...)
            end
        end # _th_readdir
        return _core_biomv, _gem_biomv
    end
    global core_biomv, gem_biomv = ret

    return
end

## --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# Correlations
let
    h0 = Histogram(
        -0.5:0.05:3.0,                  # core_biomv
        -0.5:0.05:3.0,                  # gem_biomv
    )
    @time for (core_biom, gem_biom) in zip(core_biomv, gem_biomv)
        count!(h0, (core_biom, gem_biom))
    end

    # Plot
    # 2D
    # TODO: add cor(core_biomv, gem_biomv) text
    f = Figure()
    g = GridLayout(f[1, 1])
    ax = Axis(g[1:3,2:5]; 
        title = "CORE vs GEM growth",
        limits = (-0.1, 2.5, -0.1, 2.5), 
    )
    x1 = collect(keys(h0, 1))
    x2 = collect(keys(h0, 2))
    @show length(x2)
    w = collect(values(h0))
    scatter!(ax, x1, x2; 
        colormap = :viridis, markersize = 20, 
        color = log10.(w) ./ maximum(log10, w), 
        alpha = 1.0
    )
    lines!(ax, -0.0:0.05:2.3, -0.0:0.05:2.3;
        linestyle = :dot, 
        color = :black,
        linewidth = 10
    )
    Colorbar(g[1:3, 6]; 
        label = "log10(count)",
        colormap = :viridis, limits = extrema(log10.(w)), 
    )
    
    # marginals
    h1 = marginal(h0, 1)
    ax = Axis(g[4,2:5]; 
        xlabel = "core biom [1/h]", ylabel = "count",
        limits = (-0.1, 2.5, nothing, nothing)
    )
    barplot!(ax, collect(keys(h1, 1)), collect(values(h1));
        color = :black, gap = -1
    )

    # marginals
    h1 = marginal(h0, 2)
    ax = Axis(g[1:3,1]; 
        xlabel = "count", ylabel = "gem biom [1/h]",
        limits = (nothing, nothing, -0.1, 2.5),
        xticklabelrotation = pi/4
    )
    barplot!(ax, collect(keys(h1, 1)), collect(values(h1));
        direction=:x, color = :black, gap = -1
    )

    colgap!(g, 1, -20)
    rowgap!(g, 3, -5)
    
    f
end

## --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# Comulative
let
    f =Figure()
    ax = Axis(f[1,1]; 
        title = "CORE vs GEM growth",
        xlabel = "downregulation index (sorted)",
        ylabel = "biomass [1/h]",
    )
    lines!(ax, sort(core_biomv);
        linewidth = 5, 
        label = "core"
    )
    lines!(ax, sort(gem_biomv);
        linewidth = 5,
        label = "gem"
    )
    axislegend(ax; position = :lt)
    f
end