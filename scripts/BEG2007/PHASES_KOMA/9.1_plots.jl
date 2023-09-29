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
# Collect max biomass fba solutions
let
    # context
    _simver = "ECOLI-CORE-BEG2007-PHASE_I-0.1.0"
    _load_contextdb(_simver)

    CORE_XLEP_DB = query(["ROOT", "CORE_XLEP"])
    GEM_XLEP_DB = query(["ROOT", "GEM_XLEP"])
    global RXN_MAP = query(["ROOT", "GEM_CORE_RXN_MAP"])["rxn_map"]
    
    # lep
    core_elep0 = CORE_XLEP_DB["core_elep0"][]
    global core_lep0 = lepmodel(core_elep0)
    core_elep0 = nothing
    # TDOD: boxing kill GEM growth, using unbox lep0
    global gem_lep0 = GEM_XLEP_DB["gem_lep0"][]
    global CORE_RXNI_MAP = Dict(id => i for (i, id) in enumerate(colids(core_lep0)))
    global GEM_RXNI_MAP = Dict(id => i for (i, id) in enumerate(colids(gem_lep0)))

    n0 = 0 # init file
    n1 = Inf # non-ignored file count
    cid = (@__FILE__, _simver, "biomass:core vs gem", n0, n1)
    lk = ReentrantLock()
    _, ret = withcachedat(PROJ, :get!, cid) do
        _h0 = Histogram(
            -1000.0:0.01:1000.0,                  # core_rxns
            -1000.0:0.01:1000.0,                  # core_rxns
        )
        h_th_pool = Dict()
        _th_readdir(_simver; n0, n1, nthrs = 10) do bbi, bb
            haskey(bb["meta"], "core_biomass_fba.ver") || return :ignore
            haskey(bb["meta"], "gem_biomass_fba.ver") || return :ignore
            # load frame
            feasets_db = bb["core_feasets"]
            h_pool = get!(h_th_pool, threadid(), Dict())
            _th_gem_biomv, _th_core_biomv = Float64[], Float64[]
            for feasets_blob0 in feasets_db
                for (_fealen, feasets_blob1) in feasets_blob0
                    _core_sol = feasets_blob1["core_biomass_fba.solution"]
                    _gem_sol  = feasets_blob1["gem_biomass_fba.solution"]
                    isempty(_core_sol) && continue
                    isempty(_gem_sol) && continue
                    # for (core_rxn, gem_rxn) in RXN_MAP
                    for core_rxn in colids(core_lep0)
                        gem_rxn = RXN_MAP[core_rxn]
                        core_rxni = CORE_RXNI_MAP[core_rxn]
                        gem_rxni = GEM_RXNI_MAP[gem_rxn]
                        h = get!(h_pool, core_rxn) do 
                            deepcopy(_h0)
                        end
                        count!(h, 
                            (_core_sol[core_rxni], _gem_sol[gem_rxni])
                        )
                    end
                end
            end # for feasets_blob0
        end # _th_readdir
        # reduce
        h0_pool = Dict()
        for (th, pool) in h_th_pool
            for (core_rxn, hi) in pool
                h0 = get!(h0_pool, core_rxn) do 
                    deepcopy(_h0)
                end
                merge!(h0, hi) 
            end
        end
        return h0_pool
    end
    global h0_pool = ret

    return
end

## --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# Growth Correlations Plots
let
    rxn = "BIOMASS_Ecoli_core_w_GAM"
    h0 = h0_pool[rxn]
             
    # Plot
    # 2D
    f = Figure()
    g = GridLayout(f[1, 1])
    ax = Axis(g[1:3,2:5]; 
        title = "CORE vs GEM [$(rxn)]",
        limits = (-0.1, 2.5, -0.1, 2.5), 
    )
    x1 = collect(keys(h0, 1))
    x2 = collect(keys(h0, 2))
    @show length(x2)
    @show cor(x1, x2)
    w = collect(values(h0))
    sidx = sortperm(w; rev = false)
    scatter!(ax, x1[sidx], x2[sidx]; 
        colormap = :viridis, markersize = 20, 
        color = log10.(w[sidx]) ./ maximum(log10, w), 
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
        color = :black, gap = -1, 
        width = 0.03
    )

    # marginals
    h1 = marginal(h0, 2)
    ax = Axis(g[1:3,1]; 
        xlabel = "count", ylabel = "gem biom [1/h]",
        limits = (nothing, nothing, -0.1, 2.5),
        xticklabelrotation = pi/4
    )
    barplot!(ax, collect(keys(h1, 1)), collect(values(h1));
        direction=:x, color = :black, gap = -1, 
        width = 0.04
    )

    colgap!(g, 1, -20)
    rowgap!(g, 3, -5)
    
    f
end

## --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# Correlations Coes
let
    cors = Float64[]
    bioms_cor = 0.0
    biom_id = "BIOMASS_Ecoli_core_w_GAM"
    for (rxn, h0) in h0_pool
        @show rxn
        nsamples = sum(values(h0))
        nresamples = 300_000 # 
        # @show nsamples
        scale = min(nresamples / nsamples, 1.0)
        # @show scale
        # x1 = collect(keys(h0, 1))
        x1 = resample(h0, 1; scale)
        # x2 = collect(keys(h0, 2))
        x2 = resample(h0, 2; scale)
        # @show length(x1)
        _cor = cor(x1, x2)
        isnan(_cor) && continue
        push!(cors, _cor)
        rxn == biom_id && (bioms_cor = _cor)
    end
    sort!(cors)
    
    f = Figure()
    ax = Axis(f[1,1]; 
        title = "CORE vs GEM optimal",
        xlabel = "core rxn index (sorted)", 
        ylabel = "Pearson correlation"
    )
    lines!(ax, cors; 
        linewidth = 5, 
        color = :black,
        label = "all reactions"
    )
    idx = findfirst(isequal(bioms_cor), cors)
    scatter!(ax, [idx], [bioms_cor];
        markersize = 35,
        color = :black, 
        marker = :star5, 
        label = "biom"
    )
    axislegend(ax; position = :lt)
    f
end

## --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# Biomass Comulative
let
    rxn = "BIOMASS_Ecoli_core_w_GAM"
    h0 = h0_pool[rxn]
    
    f = Figure()
    ax = Axis(f[1,1]; 
        title = "CORE vs GEM growth",
        xlabel = "downset index (sorted)",
        ylabel = "biomass [1/h]",
    )

    scale = 1e-2
    samples = resample(h0, 1; scale)
    lines!(ax, eachindex(samples) ./ scale, sort(samples);
        linewidth = 5, 
        label = "core"
    )
    samples = resample(h0, 2; scale)
    lines!(ax, eachindex(samples) ./ scale, sort(samples);
        linewidth = 5,
        label = "gem"
    )
    axislegend(ax; position = :lt)
    f
end