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
    cid = (@__FILE__, _simver, "fba:core vs gem", n0, n1)
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
    return _histogram2D_grid(h0, 1, 2;
        title = "Koma sets",
        xlabel = "core biomass", 
        ylabel = "gem biomass",
        limits = (nothing, nothing, nothing, nothing),
        dim1_bar_width = 0.03,
        dim2_bar_width = 0.05,
    )
end

## --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# Correlations Coes
let
    cors = Float64[]
    bioms_cor = 0.0
    biom_id = "BIOMASS_Ecoli_core_w_GAM"
    lk = ReentrantLock()
    @threads for (rxn, h0) in collect(h0_pool)
        @show rxn
        nsamples = sum(values(h0))
        nresamples = 30_000 # 
        # @show nsamples
        scale = min(nresamples / nsamples, 1.0)
        x1 = resample(h0, 1; scale)
        x2 = resample(h0, 2; scale)
        _cor = cor(x1, x2)
        isnan(_cor) && continue
        lock(lk) do
            push!(cors, _cor)
            rxn == biom_id && (bioms_cor = _cor)
        end
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