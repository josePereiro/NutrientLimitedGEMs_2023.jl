## ------------------------------------------------------------
@time begin
    using Random
    using MetXBase
    using MetXGEMs
    using ProjFlows
    using Statistics
    using CairoMakie
    using Base.Threads
    using Combinatorics
    using NutrientLimitedGEMs
end

# TODO: Add Graphs.jl kind of functionality for getting basic stuff
# ------------------------------------------------------------
include("1_setup_sim.jl")
include("1.1_utils.jl")

## ------------------------------------------------------------
# KOMA KO length distribution
let
    n = Inf
    cid = (@__FILE__, :LENGHT, n)
    _, h0 = withcachedat(PROJ, :set!, cid) do
        _h0 = identity_histogram(Int)
        h_pool = [deepcopy(_h0) for _ in 1:nthreads()]
        @time _foreach_obj_reg(;n) do fn, obj_reg
            @show fn
            @threads for obj in obj_reg
                count!(h_pool[threadid()], length(obj["core_koma.koset"]))
            end
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
        # yscale = Makie.pseudolog10
    )
    barplot!(ax, xs[sidxs], ws[sidxs]; label = "", color = :black)
    f
end

## ------------------------------------------------------------
# KOMA combinatorics histograms
# TODO: Add histograms to objs
# TODO: create a databatch manager (using lock files, etc)
let
    
    # build histogram
    lk = ReentrantLock()
    c = 4 # combinatiric dimension
    n = Inf # n files
    cid = (:COMB, c, n)
    _, len_h_pool = withcachedat(PROJ, :get!, cid) do
        h0 = identity_histogram(Vector{Int16})
        h_pool = Dict()
        @time _foreach_obj_reg(;n) do fn, obj_reg
            @show fn
            @threads for obj in obj_reg
                koset = obj["koset"]
                h = lock(lk) do
                    get!(h_pool, (threadid(), length(koset))) do
                        deepcopy(h0)
                    end
                end
                combs = combinations(koset, c)
                foreach(combs) do comb
                    count!(h, comb)
                end
            end
        end
        # reduce
        _len_h_pool = Dict()
        for ((th, len), h) in h_pool
            h0 = get!(_len_h_pool, len) do
                deepcopy(h0)
            end
            count!(h0, h) 
        end
        return _len_h_pool
    end

    # plot
    p = plot(; 
        title = string("comb: ", c),
        xlabel = "comb index (sorted)", 
        ylabel = "count"
    )

    # lens = 1:41
    lens = keys(len_h_pool)
    colors = colormap("Grays", maximum(lens))
    sidxs = nothing
    for l in sort(collect(lens))
        haskey(len_h_pool, l) || continue

        h0 = len_h_pool[l]
        ws = collect(counts(h0))
        @show length(ws)
        all(iszero, ws) && continue

        lock(lk) do
            # if isnothing(sidxs)
                sidxs = sortperm(ws)
                st = max(div(length(sidxs), 1000), 1)
                sidxs = sidxs[1:st:end]
            # end

            xs = range(0.0, 1.0; length = length(sidxs))
            plot!(p, ws[sidxs]; 
                label = string("len: ", l), 
                lw = 3, c = colors[l], 
                # alpha = 0.9, 
                # ylim = [0, maximum(ws)]
            )
        end
    end

    p
end

## ------------------------------------------------------------
