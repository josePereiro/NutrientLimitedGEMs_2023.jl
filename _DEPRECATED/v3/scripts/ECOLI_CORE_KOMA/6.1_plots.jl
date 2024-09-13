## ------------------------------------------------------------
@time begin
    using Random
    using MetXGEMs
    using MetXBase
    using MetXOptim
    using Statistics
    using MetXNetHub
    using CairoMakie
    using Base.Threads
    using ProgressMeter
    using Combinatorics
    using NutrientLimitedGEMs_2023
end

# ------------------------------------------------------------
include("1_setup_sim.jl")
include("1.1_utils.jl")

# ------------------------------------------------------------
# TODO: make cross script (not manual)
FVA_ALG_VER = v"0.2.0"

## ------------------------------------------------------------
# Net redundancy
let
    n = Inf
    cid = (:NET_HIST, n, FVA_ALG_VER)
    _, h0 = withcachedat(PROJ, :get!, cid) do
        _h0 = identity_histogram(UInt64)
        h_pool = [deepcopy(_h0) for _ in 1:nthreads()]
        @time _foreach_obj_reg(;n) do fn, obj_reg
            @show fn
            @threads :static for obj in obj_reg
                get(obj, "fva_ver", nothing) == FVA_ALG_VER || return
                haskey(obj, "fva.kos") || return
                # global _obj = obj
                fva_kos = sort(obj["fva.kos"]::Vector{Int64})
                count!(h_pool[threadid()], hash(fva_kos))
            end
        end
        count!(_h0, h_pool...) # reduce
        return _h0
    end

    # Plots
    xs = collect(bins(h0))
    ws = counts(h0, xs)
    @show length(ws) / sum(ws)

    sidxs = sortperm(ws)
    f = Figure()
    ax = Axis(f[1,1]; 
        xlabel = "net index (sorted)", ylabel = "count", 
        yscale = Makie.pseudolog10
    )
    lines!(ax, ws[sidxs]; label = "")
    f
end

## ------------------------------------------------------------
# BIOM vs FVA KOLEN
let
    n = Inf
    cid = (:BIOM_FVA_KOLEN, n, FVA_ALG_VER)
    _, (fvako_biom) = withcachedat(PROJ, :get!, cid) do
        fvako_biom_pool = [Dict{Tuple{UInt64, Int64}, Float64}() for _ in 1:nthreads()]
        @time _foreach_obj_reg(;n) do fn, obj_reg
            @show fn
            @threads :static for obj in obj_reg
                get(obj, "fva_ver", nothing) == FVA_ALG_VER || return
                haskey(obj, "fva.kos") || return
                haskey(obj, "biom") || return
                # global _obj = obj
                fva_kos = sort(obj["fva.kos"]::Vector{Int64})
                
                th = threadid()
                _fvako_biom = fvako_biom_pool[th]
                key = (hash(fva_kos), length(fva_kos))
                biom = obj["biom"]["feaset"]
                biom0 = get(_fvako_biom, key, -Inf)
                _fvako_biom[key] = max(biom, biom0)
            end
        end
        return merge(fvako_biom_pool...)
    end

    # Plots
    f = Figure()
    ax = Axis(f[1,1]; 
        xlabel = "fva_koset len", ylabel = "biom",
        limits = ((0, 40), (-0.1, 2)),
        # yscale = Makie.pseudolog10
    )
    len_vec = last.(keys(fvako_biom))
    biom_vec = collect(values(fvako_biom))
    valid_idxs = len_vec .!= 35
    violin!(ax, len_vec[valid_idxs], biom_vec[valid_idxs]; 
        label = "", 
        # side = :left, 
        width = 7, 
        # color = (:blue, 0.3)
    )
    f
end
