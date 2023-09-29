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
# shadown price histogram 
let
    # context
    _simver = "ECOLI-CORE-BEG2007-PHASE_I-0.1.0"
    _load_contextdb(_simver)

    n0 = 0 # init file
    n1 = Inf # non-ignored file count
    cid = (@__FILE__, _simver, "shadow price", n0, n1)
    lk = ReentrantLock()
    _, ret = withcachedat(PROJ, :set!, cid) do
        _h0 = Histogram(
            0000:10000,                        # 1 _fealen
            -1.0:0.001:1.0,                    # 2 EX_glc__D_e
            -1.0:0.001:1.0,                    # 3 EX_lac__D_e
            -1.0:0.001:1.0,                    # 4 EX_malt_e
            -1.0:0.001:1.0,                    # 5 EX_gal_e
            -1.0:0.001:1.0,                    # 6 EX_glyc_e
            -1.0:0.001:1.0,                    # 7 EX_ac_e
        )
        h_pool = Dict()
        _th_readdir(_simver; n0, n1, nthrs = 1) do bbi, bb
            haskey(bb["meta"], "core_ep.ver") || return :ignore
            haskey(bb["meta"], "core_nut_sp.ver") || return :ignore
            haskey(bb["meta"], "core_biomass_fba.ver") || return :ignore
            # load frame
            feasets_db = bb["core_feasets"]
            _h = get!(() -> deepcopy(_h0), h_pool, threadid())
            for feasets_blob0 in feasets_db
                for (_fealen, feasets_blob1) in feasets_blob0
                    # "EX_glc__D_e", "EX_lac__D_e", "EX_malt_e",
                    # "EX_gal_e", "EX_glyc_e", "EX_ac_e"
                    t = (
                        _fealen, 
                        feasets_blob1["core_nut_sp.EX_glc__D_e.obj_m"],
                        feasets_blob1["core_nut_sp.EX_lac__D_e.obj_m"],
                        feasets_blob1["core_nut_sp.EX_malt_e.obj_m"],
                        feasets_blob1["core_nut_sp.EX_gal_e.obj_m"],
                        feasets_blob1["core_nut_sp.EX_glyc_e.obj_m"],
                        # feasets_blob1["core_nut_sp.EX_ac_e.obj_m"],
                        1,
                    )
                    any(isnan, t) && continue
                    count!(_h, t)
                end
            end # for feasets_blob0
        end # _th_readdir
        # reduce
        merge!(_h0, values(h_pool)...) # reduce
        return _h0
    end
    global h0 = ret

    return
end

## --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
let
    # Plot
    # 2D
    return _histogram2D_grid(h0, 1, 4;
        title = "Koma sets",
        xlabel = "downregulation length", 
        ylabel = "shadow price",
        limits = (-10, 170, nothing, nothing),
        dim1_bar_width = 2.0,
        dim2_bar_width = 0.002,
    )
end

## --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# sp Glc > 0 only
let
    # context
    _simver = "ECOLI-CORE-BEG2007-PHASE_I-0.1.0"
    _load_contextdb(_simver)

    n0 = 0 # init file
    n1 = Inf # non-ignored file count
    cid = (@__FILE__, _simver, "shadow price", n0, n1)
    lk = ReentrantLock()
    _, ret = withcachedat(PROJ, :set!, cid) do
        _h0 = Histogram(
            0000:10000,                        # 1 _fealen
            -1.0:0.001:1.0,                    # 2 EX_glc__D_e
            -1.0:0.001:1.0,                    # 3 EX_lac__D_e
            -1.0:0.001:1.0,                    # 4 EX_malt_e
            -1.0:0.001:1.0,                    # 5 EX_gal_e
            -1.0:0.001:1.0,                    # 6 EX_glyc_e
            -1.0:0.001:1.0,                    # 7 EX_ac_e
        )
        h_pool = Dict()
        _th_readdir(_simver; n0, n1, nthrs = 1) do bbi, bb
            haskey(bb["meta"], "core_ep.ver") || return :ignore
            haskey(bb["meta"], "core_nut_sp.ver") || return :ignore
            haskey(bb["meta"], "core_biomass_fba.ver") || return :ignore
            # load frame
            feasets_db = bb["core_feasets"]
            _h = get!(() -> deepcopy(_h0), h_pool, threadid())
            for feasets_blob0 in feasets_db
                for (_fealen, feasets_blob1) in feasets_blob0
                    # "EX_glc__D_e", "EX_lac__D_e", "EX_malt_e",
                    # "EX_gal_e", "EX_glyc_e", "EX_ac_e"
                    dat = (
                        _fealen, 
                        feasets_blob1["core_nut_sp.EX_glc__D_e.obj_m"],
                        feasets_blob1["core_nut_sp.EX_lac__D_e.obj_m"],
                        feasets_blob1["core_nut_sp.EX_malt_e.obj_m"],
                        feasets_blob1["core_nut_sp.EX_gal_e.obj_m"],
                        feasets_blob1["core_nut_sp.EX_glyc_e.obj_m"],
                        # feasets_blob1["core_nut_sp.EX_ac_e.obj_m"],
                        1.0,
                    )
                    any(isnan, dat) && continue
                    abs(dat[3]) > 0.01 || continue
                    count!(_h, dat)
                end
            end # for feasets_blob0
        end # _th_readdir
        # reduce
        merge!(_h0, values(h_pool)...) # reduce
        return _h0
    end
    global h0 = ret

    return
end

## --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
let
    # Plot
    # 2D
    return _histogram2D_grid(h0, 1, 2;
        title = "Koma sets",
        xlabel = "downregulation length",
        ylabel = "shadow price",
        limits = (-10, 170, nothing, nothing),
        dim1_bar_width = 2.0,
        dim2_bar_width = 0.002,
    )
end

## --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---