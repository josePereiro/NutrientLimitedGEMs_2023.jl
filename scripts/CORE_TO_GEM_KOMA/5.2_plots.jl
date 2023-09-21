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
include("1_setup_sim.jl")
include("1.1_utils.jl")

## ------------------------------------------------------------
# identity-histogram
let
    n0 = 0
    n1 = 100
    cid = (@__FILE__, "feasets:identity-histogram", n1)
    lk = ReentrantLock()
    _, h0 = withcachedat(PROJ, :get!, cid) do
        _h0 = identity_histogram(UInt64)
        h_pool = [deepcopy(_h0) for _ in 1:nthreads()]
        _th_readdir(n1; nthrs = 10) do bbi, bb

            haskey(bb["meta"], "core_koma.ver") || return :ignore
            haskey(bb["meta"], "core_feasets.ver") || return :ignore
            haskey(bb["meta"], "core_biomass.ver") || return :ignore
            haskey(bb["meta"], "core_nut_sp.ver") || return :ignore

            # load frame
            feasets_db = bb["core_feasets"]
            core_strip_db = bb["core_strip"]

            _th_fealens, _th_feabioms, _th_feahashs = Float64[], Float64[], UInt64[]
            for (feasets_blob0, strip_blob0) in zip(feasets_db, core_strip_db)
                strip_koset = strip_blob0["koset"]
                for (_fealen, feasets_blob1) in feasets_blob0
                    _feahash = hash(sort(strip_koset[1:_fealen]))
                    count!(h_pool[threadid()], _feahash)
                end
            end
            return :continue
        end # _th_readdir
        count!(_h0, h_pool...) # reduce
        return _h0
    end
    
    # Plots
    xs = collect(bins(h0))
    ws = counts(h0, xs)
    @show length(ws) / sum(ws)

    # filter
    # sidxs = findall(w -> w > 1e2, ws)
    sidxs = eachindex(ws)
    xs = ws[sidxs]
    ws = ws[sidxs]
    # sort
    sidxs = sortperm(ws)
    f = Figure()
    ax = Axis(f[1,1]; 
        title = basename(@__DIR__),
        xlabel = "feasets.downset index", ylabel = "count",
    )
    lines!(ax, eachindex(xs[sidxs]), ws[sidxs]; label = "", color = :black)
    f
end