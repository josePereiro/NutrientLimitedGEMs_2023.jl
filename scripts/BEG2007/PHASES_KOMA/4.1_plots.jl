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
    using Base.Threads: Atomic
    using NutrientLimitedGEMs
end

# TODO: Add Graphs.jl kind of functionality for getting basic stuff
# ------------------------------------------------------------
include("1_setup.jl")
include("2_utils.jl")

## ------------------------------------------------------------
# identity histogram
let
    n1 = 10
    cid = (@__FILE__, :IDENTITY, n1)
    _, h0 = withcachedat(PROJ, :get!, cid) do
        _h0 = identity_histogram(UInt64)
        h_pool = [deepcopy(_h0) for _ in 1:nthreads()]
        
        _th_readdir(; n1, info_frec = 10) do bbi, bb
            haskey(bb["meta"], "core_koma.ver") || return false
            for blob in bb["core_koma"]
                count!(h_pool[threadid()], hash(blob["koset"]))
            end
            return false # continue
        end # _th_readdir

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
# koma.length histogram
let
    n1 = Inf
    cid = (@__FILE__, "koma.lenght histogram", n1)
    _, ret = withcachedat(PROJ, :set!, cid) do
        _h0 = Histogram(
            0:10_000,      # koma len
            0:10_000,      # rxn idx
        )
        h_pool = Dict()
        _th_readdir(; n1, info_frec = 10) do bbi, bb
            haskey(bb["meta"], "core_koma.ver") || return :ignore
            _h = get!(() -> deepcopy(_h0), h_pool, threadid())
            for blob in bb["core_koma"]
                koset = blob["koset"]
                koset_length = length(koset)
                for rxni in koset
                    count!(_h, (koset_length, rxni))
                end
            end
        end # _th_readdir
        merge!(_h0, values(h_pool)...) # reduce
        return _h0
    end
    h0 = ret

    _histogram2D_grid(h0, 1, 2;
        title = "Koma sets", 
        xlabel = "koma len", 
        ylabel = "rxn idx",
        limits = (nothing, nothing, nothing, nothing),
        dim1_bar_width = 1.6,
        dim2_bar_width = 0.6,
    )
end

## ------------------------------------------------------------
# komaset index vs rxn idx
let
    n1 = Inf
    cid = (@__FILE__, "koma.lenght vs rxn idx", n1)
    _, ret = withcachedat(PROJ, :set!, cid) do
        _h0 = Histogram(
            0:10_000,      # koma indx
            0:10_000       # rxn idx
        )
        h_pool = Dict()
        _th_readdir(; n1, info_frec = 10) do bbi, bb
            haskey(bb["meta"], "core_koma.ver") || return :ignore
            _h = get!(() -> deepcopy(_h0), h_pool, threadid())
            for blob in bb["core_koma"]
                for (i, rxni) in enumerate(blob["koset"])
                    count!(_h, (i, rxni))
                end
            end
        end # _th_readdir
        merge!(_h0, values(h_pool)...) # reduce
        return _h0
    end
    h0 = ret

    return _histogram2D_grid(h0, 1, 2;
        title = "Koma sets",
        xlabel = "set index (sorted)", 
        ylabel = "rxn index",
        limits = (nothing, nothing, nothing, nothing),
        dim1_bar_width = 0.6,
        dim2_bar_width = 0.6,
    )

end

# ## ------------------------------------------------------------
# # KOMA combinatorics histograms
# # TODO: Add histograms to objs
# # TODO: create a databatch manager (using lock files, etc)
# let
    
#     # build histogram
#     lk = ReentrantLock()
#     c = 1 # combinatiric dimension
#     n = Inf # n files
#     cid = (:COMB, c, n)
#     _, len_h_pool = withcachedat(PROJ, :get!, cid) do
#         h0 = identity_histogram(Vector{Int16})
#         h_pool = Dict()
#         @time _foreach_obj_reg(;n) do fn, obj_reg
#             @show fn
#             @threads for obj in obj_reg
#                 koset = obj["core_koma.koset"]
#                 h = lock(lk) do
#                     get!(h_pool, (threadid(), length(koset))) do
#                         deepcopy(h0)
#                     end
#                 end
#                 combs = combinations(koset, c)
#                 foreach(combs) do comb
#                     count!(h, comb)
#                 end
#             end
#         end
#         # reduce
#         _len_h_pool = Dict()
#         for ((th, len), h) in h_pool
#             h0 = get!(_len_h_pool, len) do
#                 deepcopy(h0)
#             end
#             count!(h0, h) 
#         end
#         return _len_h_pool
#     end

#     # plot
#     f = Figure()
#     ax = Axis(f[1,1]; 
#         title = string("comb: ", c),
#         xlabel = "comb index (sorted)", 
#         ylabel = "count"
#     )

#     # lens = 1:41
#     lens = len_h_pool |> keys |> collect |> sort
#     lens = lens[1:5:end]
#     # colors = colormap("Grays", maximum(lens))
#     sidxs = nothing
#     for l in sort(collect(lens))
#         haskey(len_h_pool, l) || continue

#         h0 = len_h_pool[l]
#         ws = collect(counts(h0))
#         # @show length(ws)
#         all(iszero, ws) && continue

#         lock(lk) do
#             # if isnothing(sidxs)
#                 sidxs = sortperm(ws)
#                 # st = max(div(length(sidxs), 1000), 1)
#                 # sidxs = sidxs[1:st:end]
#             # end

#             xs = range(0.0, 1.0; length = length(sidxs))
#             lines!(ax, ws[sidxs]; 
#                 label = string("len: ", l), 
#                 # lw = 3, 
#                 # color = colors[l], 
#                 # alpha = 0.9, 
#                 # ylim = [0, maximum(ws)]
#             )
#         end
#     end

#     f
# end

# ## ------------------------------------------------------------


