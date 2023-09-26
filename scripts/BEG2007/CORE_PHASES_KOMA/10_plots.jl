## ------------------------------------------------------------
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

# TODO: Add Graphs.jl kind of functionality for getting basic stuff
# ------------------------------------------------------------
include("1_setup.jl")
include("2_utils.jl")

## ------------------------------------------------------------
# biomass vs downset length histogram
let
    n0 = 0 # init file
    n1 = 3 # non-ignored file count
    cid = (@__FILE__, "feasets:length+biomass+count", n0, n1)
    lk = ReentrantLock()
    _, h0 = withcachedat(PROJ, :get!, cid) do
        _h0 = Histogram(
            0:5000,      # feasets len
            0.0:0.05:3.0, # biom
        )
        h_pool = [deepcopy(_h0) for _ in 1:nthreads()]
        _th_readdir(n1, n0; nthrs = 1) do bbi, bb

            haskey(bb["meta"], "core_koma.ver") || return :ignore
            haskey(bb["meta"], "core_feasets.ver") || return :ignore
            haskey(bb["meta"], "core_biomass.ver") || return :ignore
            haskey(bb["meta"], "core_nut_sp.ver") || return :ignore

            # load frame
            feasets_db = bb["core_feasets"]

            for feasets_blob0 in feasets_db
                for (_fealen, feasets_blob1) in feasets_blob0
                    _biom = feasets_blob1["core_biomass.biom"]
                    count!(h_pool[threadid()], (_fealen, _biom))
                end
            end
            return :continue
        end # _th_readdir
        merge!(_h0, h_pool...) # reduce
        return _h0
    end

    # Plot
    # 2D
    f = Figure()
    g = GridLayout(f[1, 1])
    ax = Axis(g[1:3,2:5]; 
        title = basename(@__DIR__),
        limits = (-10, 200, -0.2, 2.6)
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
    Colorbar(g[1:3, 6]; 
        label = "log10(count)",
        colormap = :viridis, limits = extrema(log10.(w)), 
    )
    
    # marginals
    h1 = marginal(h0, 1)
    ax = Axis(g[4,2:5]; 
        xlabel = "downset length", ylabel = "count",
        limits = (-10, 200, nothing, nothing, )
    )
    barplot!(ax, collect(keys(h1, 1)), collect(values(h1));
        color = :black, gap = -1
    )

    # marginals
    h1 = marginal(h0, 2)
    ax = Axis(g[1:3,1]; 
        xlabel = "count", ylabel = "biom",
        limits = (nothing, nothing, -0.2, 2.6),
        xticklabelrotation = pi/4
    )
    barplot!(ax, collect(keys(h1, 1)), collect(values(h1));
        direction=:x, color = :black, gap = -1
    )

    colgap!(g, 1, -20)
    rowgap!(g, 3, -5)
    
    f
end

## ------------------------------------------------------------
# feaset lenght vs shadow price of nutrient
let
    n0 = 0 # init file
    n1 = Inf # non-ignored file count
    exch = ["EX_glc__D_e", "EX_lac__D_e", "EX_malt_e",
    "EX_gal_e", "EX_glyc_e", "EX_ac_e", "EX_ac_e"][7]
    lk = ReentrantLock()
    cid = (@__FILE__, "feasets:nutrient.price vs downset.lenght", n0, n1, exch)
    _, h0 = withcachedat(PROJ, :get!, cid) do
        _h0 = Histogram(
            0:5000,          # feasets len
            -10.0:0.0025:10.0, # shadow price
        )
        h_pool = [deepcopy(_h0) for _ in 1:nthreads()]
        _th_readdir(n1, n0; nthrs = 10) do bbi, bb

            # filter
            haskey(bb["meta"], "core_koma.ver") || return :ignore
            haskey(bb["meta"], "core_feasets.ver") || return :ignore
            haskey(bb["meta"], "core_biomass.ver") || return :ignore
            haskey(bb["meta"], "core_nut_sp.ver") || return :ignore

            # load frame
            feasets_db = bb["core_feasets"]

            # _th_fealens, _th_obj_ms = Float64[], Float64[]
            for feasets_blob0 in feasets_db
                for (_fealen, feasets_blob1) in feasets_blob0
                    mkey = string("core_nut_sp.", exch, ".obj_m")
                    _price = feasets_blob1[mkey]
                    # abs(_price) < 1e-3 && continue # ignore non exch limited
                    count!(h_pool[threadid()], (_fealen, _price))
                end
            end

            return :continue
        end # _th_readdir
        merge!(_h0, h_pool...) # reduce
        return _h0
    end

    # Plot
    # 2D
    f = Figure()
    g = GridLayout(f[1, 1])
    ax = Axis(g[1:3,2:5]; 
        title = string(
            basename(@__DIR__), "\n",
            "exch: ", exch
        ),
        # limits = (-10, 200, -0.2, 2.6)
    )
    x1 = collect(keys(h0, 1)) # feaslen
    x2 = collect(keys(h0, 2)) # price
    w = collect(values(h0))
    #  TODO: filter histogram
    # fidxs = findall((x) -> abs(x) > 1e-2, x2)
    # x1 = x1[fidxs]
    # x2 = x2[fidxs]
    # w = w[fidxs]
    scatter!(ax, x1, x2; 
        colormap = :viridis, markersize = 20, 
        color = log10.(w) ./ maximum(log10, w), 
        alpha = 1.0
    )
    Colorbar(g[1:3, 6]; 
        label = "log10(count)",
        colormap = :viridis, limits = extrema(log10.(w)), 
    )
    
    # marginals
    h1 = marginal(h0, 1)
    ax = Axis(g[4,2:5]; 
        xlabel = "downset length", ylabel = "count",
        # limits = (-10, 200, nothing, nothing, )
    )
    barplot!(ax, collect(keys(h1, 1)), collect(values(h1));
        color = :black, gap = -1
    )

    # marginals
    h1 = marginal(h0, 2)
    ax = Axis(g[1:3,1]; 
        xlabel = "count", ylabel = "shadow price",
        # limits = (nothing, nothing, -0.2, 2.6),
        xticklabelrotation = pi/4
    )
    barplot!(ax, collect(keys(h1, 1)), collect(values(h1));
        direction=:x, color = :black, gap = -1
    )

    colgap!(g, 1, -15)
    rowgap!(g, 3, -5)
    
    f
end

## ------------------------------------------------------------
## ------------------------------------------------------------
