## ------------------------------------------------------------
@time begin
    using CairoMakie
    using Random
    using MetXGEMs
    using MetXBase
    using MetXOptim
    using Statistics
    using MetXNetHub
    using Base.Threads
    using ProgressMeter
    using Combinatorics
    using NutrientLimitedGEMs
end

# TODO: Add Graphs.jl kind of functionality for getting basic stuff
# ------------------------------------------------------------
include("1_setup_sim.jl")
include("1.1_utils.jl")

## ------------------------------------------------------------
# Grow Histograms
let
    n = Inf
    nut_id = "EX_lac__D_e"
    cid = (:BIOM, n, nut_id)
    _, h0 = withcachedat(PROJ, :get!, cid) do
        _h0 = identity_histogram(Float64)
        h_pool = [deepcopy(_h0) for _ in 1:nthreads()]
        # nut_ids = ["EX_glc__D_e", "EX_lac__D_e", "EX_ac_e"]
        @time _foreach_obj_reg(;n) do fn, obj_reg
            @show fn
            @threads for obj in obj_reg
                haskey(obj, "biom") || continue
                biom_reg = obj["biom"]
                # for nut_id in nut_ids
                    count!(h_pool[threadid()], biom_reg[nut_id])
                # end
            end
        end
        count!(_h0, h_pool...) # reduce
        return _h0
    end

    # plot
    xs = collect(x for x in bins(h0) if x > 0.01)
    ws = counts(h0, xs)
    sidxs = sortperm(xs)
    p = plot(;title = nut_id)
    bar!(p, xs[sidxs], ws[sidxs]; label = "")
end

## ------------------------------------------------------------

## ------------------------------------------------------------
# Nut essential length histograms
let
    n = Inf
    ps = Plots.Plot[]
    batch_size = 3
    nut_ids = ["EX_glc__D_e", "EX_lac__D_e", "EX_ac_e"]
    p = plot(; title = "ecoli core", xlabel = "feasible set length", ylabel = "count")
    for nut_id in nut_ids
        cid = (:ESS_LENGHT, n, nut_id)
        _, h0 = withcachedat(PROJ, :get!, cid) do
            _h0 = identity_histogram(Int)
            @time _foreach_obj_reg(;n) do fn, obj_reg
                @show fn
                for obj in obj_reg
                    haskey(obj, "biom") || continue
                    haskey(obj, "koset") || continue
                    koset = obj["koset"]
                    length(koset) == batch_size && continue
                    biom_reg = obj["biom"]
                    biom_reg[nut_id] < 0.1 || continue
                    count!(_h0, length(koset) - batch_size)
                end
            end
            return _h0
        end
        xs = collect(bins(h0))
        ws = counts(h0, xs)
        sidxs = sortperm(xs)
        plot!(xs[sidxs], ws[sidxs]; label = nut_id)
    end
    p
    # plot(ps...)
end

## ------------------------------------------------------------
# BIOM vs KOLEN
let
    n = Inf
    cid = (:BIOM_KOLEN, n)
    batch_size = 3
    _, (len_vec, biom_vec) = withcachedat(PROJ, :get!, cid) do
        len_vec_pool = [Int64[] for _ in 1:nthreads()]
        biom_vec_pool = [Float64[] for _ in 1:nthreads()]
        @time _foreach_obj_reg(;n) do fn, obj_reg
            @show fn
            @threads for obj in obj_reg
                haskey(obj, "biom") || return
                # global _obj = obj
                feaset = obj["koset"][1:(end-batch_size)]
                biom = obj["biom"]["feaset"]
                th = threadid()
                push!(len_vec_pool[th], length(feaset))
                push!(biom_vec_pool[th], biom)
            end
        end
        return vcat(len_vec_pool...), vcat(biom_vec_pool...)
    end

    # Plots
    f = Figure()
    ax = Axis(f[1,1]; 
        xlabel = "koset len", ylabel = "biom",
        limits = ((0, 34), (-0.1, 2)),
    )
    violin!(ax, len_vec, biom_vec; 
        label = "", 
        width = 5, 
        # side = :left, 
    )
    f
end
