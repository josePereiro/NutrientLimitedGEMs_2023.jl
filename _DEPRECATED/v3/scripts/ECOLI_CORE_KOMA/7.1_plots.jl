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
    using NutrientLimitedGEMs
end

# ------------------------------------------------------------
include("1_setup_sim.jl")
include("1.1_utils.jl")

# ------------------------------------------------------------
# TODO: make cross script (not manual)
STRIP_ALG_VER = v"0.1.0"

## ------------------------------------------------------------
# Striped lenght histogram
let
    n = Inf
    cid = (@__FILE__, :NET_HIST, n, STRIP_ALG_VER)
    _, h0 = withcachedat(PROJ, :get!, cid) do
        _h0 = identity_histogram(Int64)
        _unique_set = Set{UInt64}()
        lk = ReentrantLock()
        h_pool = [deepcopy(_h0) for _ in 1:nthreads()]
        @time _foreach_obj_reg(;n) do fn, obj_reg
            @show fn
            @threads :static for obj in obj_reg
                get(obj, "strip_ver", nothing) == STRIP_ALG_VER || continue
                haskey(obj, "strip.koset") || continue
                # global _obj = obj
                strip_koset = obj["strip.koset"]
                koset_hash = hash(strip_koset)
                koset_hash ∈ _unique_set && continue
                lock(lk) do
                    push!(_unique_set, koset_hash)
                end
                count!(h_pool[threadid()], length(strip_koset))
            end
            @show length(_unique_set)
        end
        count!(_h0, h_pool...) # reduce
        return _h0
    end

    # Plots
    xs = collect(bins(h0))
    ws = counts(h0, xs)
    @show length(ws) / sum(ws)

    sidxs = sortperm(xs)
    f = Figure()
    ax = Axis(f[1,1]; 
        title = basename(@__DIR__),
        xlabel = "strip.koset length", ylabel = "count", 
        # yscale = Makie.pseudolog10
    )
    barplot!(ax, xs[sidxs], ws[sidxs]; label = "", color = :black)
    f
end

## ------------------------------------------------------------
# Biom vs unique strip.koset
let
    biom_id = "feaset"
    n = 5
    cid = (@__FILE__, biom_id, :BIOM_HIST, n, STRIP_ALG_VER)
    _, h0 = withcachedat(PROJ, :get!, cid) do
        _h0 = identity_histogram(Float64)
        _unique_set = Set{UInt64}()
        lk = ReentrantLock()
        h_pool = [deepcopy(_h0) for _ in 1:nthreads()]
        @time _foreach_obj_reg(;n) do fn, obj_reg
            @show fn
            @threads :static for obj in obj_reg
                get(obj, "strip_ver", nothing) == STRIP_ALG_VER || continue
                haskey(obj, "strip.koset") || continue
                haskey(obj, "biom") || continue
                # global _obj = obj
                strip_koset = obj["strip.koset"]
                koset_hash = hash(strip_koset)
                koset_hash ∈ _unique_set && continue
                lock(lk) do
                    push!(_unique_set, koset_hash)
                end
                count!(h_pool[threadid()], obj["biom"][biom_id])
            end
            @show length(_unique_set)
        end # _foreach_obj_reg
        count!(_h0, h_pool...) # reduce
        return _h0
    end # cache

    # Plots
    xs = collect(bins(h0))
    ws = counts(h0, xs)
    @show length(ws) / sum(ws)

    sidxs = sortperm(xs)
    f = Figure()
    ax = Axis(f[1,1]; 
        title = basename(@__DIR__),
        xlabel = "strip.koset length", ylabel = string("biom (", biom_id, ")"), 
        # yscale = Makie.pseudolog10
    )
    barplot!(ax, xs[sidxs], ws[sidxs]; label = "", color = :black)
    f
end