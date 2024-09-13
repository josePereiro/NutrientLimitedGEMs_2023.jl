## --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
@time begin
    using MetXGEMs
    using MetXBase
    using Statistics
    using ProjFlows
    using CairoMakie
    using BlobBatches
    using Base.Threads
    using NutrientLimitedGEMs_2023
end

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
include("1_setup.jl")
include("2_utils.jl")

## --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# Fva range/center
## --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
let
    # context
    _simver = "ECOLI-CORE-BEG2007-PHASE_0"
    _load_contextdb(_simver)

    n0 = 0 # init file
    n1 = 100 # non-ignored file count
    cid = (@__FILE__, _simver, "fva: range/center moments", n0, n1)
    lk = ReentrantLock()
    _, ret = withcachedat(PROJ, :set!, cid) do
        _h0 = Histogram(
            0:1:10000,                        #= _fealen =#
            0e0:1e0:1e9,                        #= mean =#
            0e0:1e0:1e9,                        #= std =#
            0e0:1e0:1e9,                        #= mean =#
            0e0:1e0:1e9,                        #= std =#
        )
        h_pool = Dict()
        _th_readdir(_simver; n0, n1, nthrs = 10) do bbi, bb
            haskey(bb["meta"], "core_fva.ver") || return :ignore
            # load frame
            feasets_db = bb["core_feasets"]
            _h = get!(() -> deepcopy(_h0), h_pool, threadid())
            for feasets_blob0 in feasets_db
                for (_fealen, feasets_blob1) in feasets_blob0
                    _fvalb = Float64.(feasets_blob1["core_fva.fvalb"])
                    _fvaub = Float64.(feasets_blob1["core_fva.fvaub"])
                    _range = _fvaub .- _fvalb
                    _center = _fvaub .+ _fvalb
                    range_mean = mean(_range)
                    range_std = std(_range)
                    center_mean = mean(_center)
                    center_std = std(_center)
                    count!(_h, (
                        _fealen, 
                        range_mean, range_std, 
                        center_mean, center_std
                    ))
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
        title = "Fva",
        xlabel = "downregulation set length",
        ylabel = "rxn range mean",
        # limits = (-10, 170, nothing, nothing),
        dim1_bar_width = 2.0,
        dim2_bar_width = 3.0,
    )
end

## --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
let
    # Plot
    # 2D
    return _histogram2D_grid(h0, 2, 3;
        title = "Fva",
        xlabel = "rxn range mean",
        ylabel = "rxn range std",
        # limits = (-10, 170, nothing, nothing),
        dim1_bar_width = 2.0,
        dim2_bar_width = 3.0,
    )
end

## --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
let
    # Plot
    # 2D
    return _histogram2D_grid(h0, 4, 5;
        title = "Fva",
        xlabel = "rxn center mean",
        ylabel = "rxn center std",
        # limits = (-10, 170, nothing, nothing),
        dim1_bar_width = 2.0,
        dim2_bar_width = 3.0,
    )
end

## --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# Fva range/center per rxn
## --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
let
    # context
    _simver = "ECOLI-CORE-BEG2007-PHASE_1"
    _load_contextdb(_simver)

    n0 = 0 # init file
    n1 = Inf # non-ignored file count
    cid = (@__FILE__, _simver, "fva: range/center per rxn", n0, n1)
    lk = ReentrantLock()
    _, ret = withcachedat(PROJ, :set!, cid) do
        _h0 = Histogram(
            0:1:10000,                              #= rxni =#
            -1e2:5e-1:1e2,                          #= _range =#
            -1e2:5e-1:1e2,                          #= _center =#
        )
        h_pool = Dict()
        _th_readdir(_simver; n0, n1, nthrs = 1) do bbi, bb
            haskey(bb["meta"], "core_fva.ver") || return :ignore
            # load frame
            feasets_db = bb["core_feasets"]
            _h = get!(() -> deepcopy(_h0), h_pool, threadid())
            for feasets_blob0 in feasets_db
                for (_fealen, feasets_blob1) in feasets_blob0
                    _fvalb = Float64.(feasets_blob1["core_fva.fvalb"])
                    _fvaub = Float64.(feasets_blob1["core_fva.fvaub"])
                    _range = _fvaub .- _fvalb
                    _center = _fvaub .+ _fvalb
                    for rxni in eachindex(_range)
                        count!(_h, (
                            rxni, 
                            _range[rxni], _center[rxni], 
                        ))
                    end
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
        title = "Fva",
        xlabel = "rxn index",
        ylabel = "rxn range",
        # limits = (-10, 170, nothing, nothing),
        dim1_bar_width = 1.0,
        dim2_bar_width = 1.5,
    )
end

## --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
let
    # Plot
    # 2D
    return _histogram2D_grid(h0, 1, 3;
        title = "Fva",
        xlabel = "rxn index",
        ylabel = "rxn center",
        # limits = (-10, 170, nothing, nothing),
        dim1_bar_width = 1.0,
        dim2_bar_width = 1.5,
    )
end

## --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# Fva dist from net0
## --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
let
    # context
    _simver = "ECOLI-CORE-BEG2007-PHASE_1"
    _load_contextdb(_simver)

    # lep
    CORE_XLEP_DB = query(["ROOT", "CORE_XLEP"])
    core_elep0 = CORE_XLEP_DB["core_elep0"][]
    _lb0, _ub0 = lb(core_elep0), ub(core_elep0)

    n0 = 0 # init file
    n1 = Inf # non-ignored file count
    cid = (@__FILE__, _simver, "fva: square error", n0, n1)
    lk = ReentrantLock()
    _, ret = withcachedat(PROJ, :set!, cid) do
        _h0 = Histogram(
            0:1:10000,                              #= _fealen =#
            -1e9:1e1:1e9,                          #= _mean_err =#
            -1e9:1e1:1e9,                          #= _std_err =#
        )
        h_pool = Dict()
        _th_readdir(_simver; n0, n1, nthrs = 1) do bbi, bb
            haskey(bb["meta"], "core_fva.ver") || return :ignore
            # load frame
            feasets_db = bb["core_feasets"]
            _h = get!(() -> deepcopy(_h0), h_pool, threadid())
            for feasets_blob0 in feasets_db
                for (_fealen, feasets_blob1) in feasets_blob0
                    _fvalb = Float64.(feasets_blob1["core_fva.fvalb"])
                    _fvaub = Float64.(feasets_blob1["core_fva.fvaub"])
                    _qerr = (_lb0 .- _fvalb).^2 .+ (_ub0 .- _fvaub).^2
                    _err = mean(_qerr)
                    _std = std(_qerr)
                    count!(_h, (
                        _fealen, 
                        _err, _std
                    ))
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
    return _histogram2D_grid(h0, 2, 3;
        title = "Fva",
        xlabel = "log mean sqerr",
        ylabel = "log std sqerr",
        dim1_T = v -> log10.(v),
        dim2_T = v -> log10.(v),
        limits = (2.5, 5.3, 3.5, 5.5),
        dim1_bar_width = 0.05,
        dim2_bar_width = 0.05,
        colgap = -5,
        rowgap = -5
    )
end

## --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# let
#     # Plot
#     # 2D
#     return _histogram2D_grid(h0, 1, 6;
#         title = "Fva",
#         xlabel = "rxn index",
#         ylabel = "rxn center",
#         # limits = (-10, 170, nothing, nothing),
#         dim1_bar_width = 1.0,
#         dim2_bar_width = 1.0,
#     )
# end