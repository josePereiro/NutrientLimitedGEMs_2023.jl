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
    using NutrientLimitedGEMs
end

# TODO: Add Graphs.jl kind of functionality for getting basic stuff
# ------------------------------------------------------------
include("1_setup_sim.jl")
include("1.1_utils.jl")

## ------------------------------------------------------------
# identity histogram
let
    n = Inf
    cid = (@__FILE__, :IDENTITY, n)
    _, h0 = withcachedat(PROJ, :set!, cid) do
        _h0 = identity_histogram(UInt64)
        h_pool = [deepcopy(_h0) for _ in 1:nthreads()]
        batches = readdir(BlobBatch, procdir(PROJ, [SIMVER]))
        for (bbi, bb) in enumerate(batches)
            haskey(bb["meta"], "core_koma.ver") || continue
            @show basename(rootdir(bb))
            @threads for blob in bb["core_koma"]
                count!(h_pool[threadid()], hash(blob["koset"]))
            end
            bbi > n && break
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
        title = basename(@__DIR__),
        xlabel = "koset index", ylabel = "count",
    )
    lines!(ax, eachindex(xs[sidxs]), ws[sidxs]; label = "", color = :black)
    f
end


## ------------------------------------------------------------
# koma.lenght histogram
let
    n = Inf
    cid = (@__FILE__, :LENGHT, n)
    _, h0 = withcachedat(PROJ, :set!, cid) do
        _h0 = identity_histogram(Int)
        h_pool = [deepcopy(_h0) for _ in 1:nthreads()]
        batches = readdir(BlobBatch, procdir(PROJ, [SIMVER]))
        for (bbi, bb) in enumerate(batches)
            haskey(bb["meta"], "core_koma.ver") || continue
            @show basename(rootdir(bb))
            @threads for blob in bb["core_koma"]
                count!(h_pool[threadid()], length(blob["koset"]))
            end
            bbi > n && break
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
        xlabel = "core_koma.koset length", ylabel = "count",
    )
    lines!(ax, xs[sidxs], ws[sidxs]; label = "", color = :black)
    f
end

## ------------------------------------------------------------
## ------------------------------------------------------------

## ------------------------------------------------------------
