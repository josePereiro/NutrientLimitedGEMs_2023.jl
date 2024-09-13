## ------------------------------------------------------------
@time begin
    using Random
    using MetXBase
    using MetXGEMs
    using ProjFlows
    using Statistics
    using CairoMakie
    using BlobBatches
    using Distributions
    using Base.Threads
    using Combinatorics
    using Base.Threads: Atomic
    using NutrientLimitedGEMs_2023
end

# ------------------------------------------------------------
include("1_setup.jl")
include("2_utils.jl")

## ------------------------------------------------------------
# General filter
let
    # context
    _simver = "ECOLI-CORE-BEG2007-PHASE_0"
    _load_contextdb(_simver)

    n0 = 0 # init file
    n1 = 200 # non-ignored file count
    cid = (@__FILE__, "general filter", n0, n1)
    _, ret = withcachedat(PROJ, :set!, cid) do
        _h0 = Histogram(
                   0:10000,        # _fealen
            -10000.0:0.1:10000.0,  # entropy
            -10000.0:0.1:10000.0,  # free energy
               -10.0:0.01:10.0,    # biomass
            ;
            names = ["fealen", "entropy", "free_energy", "biomass"]
        )
        h_pool = Dict()
        _th_readdir(_simver; n0, n1, nthrs = 10) do bbi, bb

            # meta filters
            haskey(bb["meta"], "core_ep.ver") || return :ignore
            haskey(bb["meta"], "core_nut_sp.ver") || return :ignore
            haskey(bb["meta"], "core_biomass_fba.ver") || return :ignore
            
            # load frame
            feasets_db = bb["core_feasets"]
            _h = get!(() -> deepcopy(_h0), h_pool, threadid())
            for feasets_blob0 in feasets_db
                for (_fealen, feasets_blob1) in feasets_blob0
                    # feaset filters
                    feasets_blob1["core_ep.status"] == :error && continue
                    _H = feasets_blob1["core_ep.entropy"]
                    _F = feasets_blob1["core_ep.free_energy"]
                    _z = feasets_blob1["core_biomass_fba.biom"]
                    glc_sp = feasets_blob1["core_nut_sp.EX_glc__D_e.obj_m"]

                    # value filters
                    any(isnan, [_H, _F, _z]) && continue
                    abs(_F) < 1000 || continue
                    
                    _z > 0.8 || continue
                    _z < 1.3 || continue

                    abs(glc_sp) > 1e-6 || continue

                    # _F > 60 || continue

                    count!(_h, (_fealen, _H, _F, _z))
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

## ------------------------------------------------------- 
# 
let
    # Plot
    # 2D
    return _histogram2D_grid(h0, 3, 4;
        xlabel = "free energy", 
        # ylabel = "growth",
        limits = (nothing, nothing, nothing, nothing),
        dim1_bar_width = 1.0,
        dim2_bar_width = 0.01,
    )
end  