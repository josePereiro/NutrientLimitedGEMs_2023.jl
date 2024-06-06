## ------------------------------------------------------------
@time begin
    using Plots
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

# ------------------------------------------------------------
include("1_setup_sim.jl")
include("1.1_utils.jl")

## ------------------------------------------------------------
# Tools
function _load_koma_regs(n = Inf)
    reg = Dict{Vector{Int16}, Symbol}()
    regfns = filter(readdir(procdir(PROJ, [SIMVER]); join = true)) do fn
        contains(basename(fn), "koma_reg")
    end
    for (i, fn) in enumerate(regfns)
        i > n && break
        regi = deserialize(fn)[:dat]
        merge!(reg, regi)
    end
    return reg
end

function _load_kosets(n = Inf)
    reg = _load_koma_regs(n)
    kosets = collect(keys(reg))
    foreach(sort!, kosets)
    unique!(kosets)
    return kosets
end

## ------------------------------------------------------------
# KOMA KO length distribution
let
    # sort!
    @time global kosets = _load_kosets()

    h = identity_histogram(Int)
    foreach(kosets) do koset
        count!(h, length(koset))
    end

    # plot
    xs = sort(collect(keys(h)))
    ws = counts(h, xs)

    p = plot(; xlabel = "koset length", ylabel = "count")
    bar!(p, xs, ws; label = "", c = :black)

end

## ------------------------------------------------------------
# KOMA combinatorics histograms
let
    # # xlep
    # xlep_db = query(["ROOT", "XLEP"])
    # global elep0 = xlep_db["elep0"][]
    # global lep0 = lepmodel(elep0)

    # sort!
    @time global kosets = _load_kosets()
    
    # loop params
    n = 3
    lens = 5:10:50
    
    # loop tools
    colors = colormap("Grays", maximum(lens))
    sidxs = nothing
    lk = ReentrantLock()

    p = plot(; 
        title = string("comb: ", n),
        xlabel = "comb index (sorted)", 
        ylabel = "count"
    )

    for l in collect(lens)

        # build histogram
        h0 = identity_histogram(Vector{Int16})
        h_pool = [deepcopy(h0) for th in 1:4]
        @threads :static for koset in kosets
            length(koset) == l || continue
            
            h = h_pool[threadid()]
            combs = combinations(koset, n)
            foreach(combs) do comb
                count!(h, comb)
            end
        end
        count!(h0, h_pool...) # reduce

        # plot
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
## ------------------------------------------------------------
## ------------------------------------------------------------
## ------------------------------------------------------------
## ------------------------------------------------------------
## ------------------------------------------------------------
# KOMA n=1 histograms
let
    # xlep
    xlep_db = query(["ROOT", "XLEP"])
    elep0 = xlep_db["elep0"][]
    lep0 = lepmodel(elep0)

    # sort!
    @time global kosets = _load_kosets()

    p = plot(; 
        xlabel = "rxn index (sorted)", 
        ylabel = "count"
    )
    lens = 50:100
    colors = colormap("Grays", maximum(lens))
    sidxs = nothing
    for l in lens
        xs = elep0.idxi
        ws = zeros(length(xs))
        for koset in kosets
            length(koset) == l || continue
            _histogram!(xs, ws, koset)
        end
        all(iszero, ws) && continue
        if isnothing(sidxs)
            sidxs = sortperm(ws)
        end
        plot!(p, ws[sidxs]; 
            label = l, 
            lw = 3, c = colors[l], 
            alpha = 0.9
        )
    end
    p

end