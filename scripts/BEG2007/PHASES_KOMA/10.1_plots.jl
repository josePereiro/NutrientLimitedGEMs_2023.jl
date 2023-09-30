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
    _simver = "ECOLI-CORE-BEG2007-PHASE_1"
    _load_contextdb(_simver)

    n0 = 0 # init file
    n1 = Inf # non-ignored file count
    cid = (@__FILE__, _simver, "ep:entropy + free energy", n0, n1)
    lk = ReentrantLock()
    _, ret = withcachedat(PROJ, :set!, cid) do
        _h0 = Histogram(
                   0:10000,                        # _fealen
            -10000.0:0.1:10000.0,                  # fealen
            -10000.0:0.1:10000.0,                  # free energy
               -10.0:0.01:10.0,                    # biomass
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
                    feasets_blob1["core_ep.status"] == :error && continue
                    _H = feasets_blob1["core_ep.entropy"]
                    _F = feasets_blob1["core_ep.free_energy"]
                    _z = feasets_blob1["core_biomass_fba.biom"]
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

## --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# Entropy vs free Enery
let
    # Plot
    # 2D
    return _histogram2D_grid(h0, 2, 3;
        title = "Entropy vs Free Energy",
        xlabel = "entropy", 
        ylabel = "free energy",
        limits = (nothing, nothing, nothing, nothing),
        dim1_bar_width = 3.0,
        dim2_bar_width = 1.0,
    )
end  

## --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# Free Energy vs len
let

    XLEP_DB = query(["ROOT", "CORE_EP.LEP0"])
    H0 = XLEP_DB["entropy"]
    F0 = XLEP_DB["free_energy"]

    return _histogram2D_grid(h0, 1, 3;
        title = "Space volume",
        xlabel = "downset lenght",
        ylabel = "~ log Vol/Vol0",
        limits = (nothing, nothing, nothing, nothing),
        dim2_T = v -> F0 .- v,
        dim1_bar_width = 2.0,
        dim2_bar_width = 1.0,
    )
end    

## --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# Free Energy vs Biomass
let

    XLEP_DB = query(["ROOT", "CORE_EP.LEP0"])
    H0 = XLEP_DB["entropy"]
    F0 = XLEP_DB["free_energy"]

    # Plot
    # 2D
    return _histogram2D_grid(h0, 4, 3;
        title = "Space volume",
        xlabel = "max biomass",
        ylabel = "~ log Vol/Vol0",
        limits = (nothing, nothing, nothing, nothing),
        dim2_T = v -> F0 .- v,
        dim1_bar_width = 0.03,
        dim2_bar_width = 4.0,
    )

end    