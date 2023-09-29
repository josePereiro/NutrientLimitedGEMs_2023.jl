## ------------------------------------------------------------
function _merge_metinfo!(dest_net1, id1, src_net2, id2)
    meti1 = metindex(dest_net1, id1)
    meti2 = metindex(src_net2, id2)
    
    for f in [:metNames, :metFormulas]
        val2 = _getindex_or_nothing(getfield(src_net2, f), meti2)
        isnothing(val2) || _setindex_or_nothing!(
            getfield(dest_net1, f), meti1, val2
        )
    end
end

function _merge_rxninfo!(dest_net1, id1, src_net2, id2)
    rxni1 = rxnindex(dest_net1, id1)
    rxni2 = rxnindex(src_net2, id2)
    
    for f in [:subSystems, :rxnNames]
        val2 = _getindex_or_nothing(getfield(src_net2, f), rxni2)
        isnothing(val2) || _setindex_or_nothing!(
            getfield(dest_net1, f), rxni1, val2
        )
    end
end

# ------------------------------------------------------------
function _sync_koma_hashs!(koma_hashs)
    # up state
    lock(PROJ) do
        # koma_hashs
        fn = procdir(PROJ, [SIMVER], "koma_hashs.jls")
        _, _koma_hashs = ldat(fn) do 
            UInt64[]
        end
        push!(koma_hashs, setdiff(_koma_hashs, koma_hashs)...)
        unique!(koma_hashs)
        sort!(koma_hashs)
        sdat(koma_hashs, fn)

        # write blobs
        println("[", getpid(), ".", threadid(), "] ", "SYNC KOMA_HASH")
    end
    return koma_hashs
end

## ------------------------------------------------------------
function _with_kos(f::Function, model, kos::Vector; zero = 0.0)
    kos = colindex(model, kos)
    lb0, ub0 = bounds(model, kos)
    try; bounds!(model, kos, zero, zero); f()
        finally; bounds!(model, kos, lb0, ub0)
    end
end

function _apply_downreg!(model, todown, DOWNREG_FACTOR)
    l, u = bounds(model, todown)
    bounds!(model, todown, l * DOWNREG_FACTOR, u * DOWNREG_FACTOR)
end
function _apply_downreg!(model, todownv::Vector, DOWNREG_FACTOR) 
    foreach(todownv) do todown    
        _apply_downreg!(model, todown, DOWNREG_FACTOR)
    end
end

function _with_downreg(f::Function, model, todownv::Vector, DOWNREG_FACTOR)
    todownv = colindex(model, todownv)
    lb0, ub0 = bounds(model, todownv)
    try; _apply_downreg!(model, todownv, DOWNREG_FACTOR); f()
    finally; bounds!(model, todownv, lb0, ub0)
    end
end

## ------------------------------------------------------------
import Base.Threads: Atomic
function _th_readdir(f::Function, sim; 
        n1 = Inf, n0 = 1,
        info_frec = 1.0, nthrs = 10, verbose = true
    )
    batches = readdir(BlobBatch, procdir(PROJ, [sim]))
    nread = Atomic{Int}(0)
    bbi = Atomic{Int}(0)
    t0 = Atomic{Float64}(-1.0)
    @threads for _ in 1:nthrs
        for bb in batches
            nread[] += 1
            nread[] >= n0 || continue

            flag = f(bbi[], bb)
            flag === :break && break
            flag === :ignore && continue

            # info
            if verbose && time() > t0[]
                println(
                    "[", getpid(), ".", threadid(), "]", 
                    "\tbbi: ", bbi[], ", nread: ", nread[]
                )
                t0[] = time() + info_frec
            end
            
            bbi[] += 1
            bbi[] > n1 && break
        end
    end
    return nothing
end
_th_readdir(f::Function; kwargs...) = _th_readdir(f::Function, SIMVER; kwargs...)

# ------------------------------------------------------------
_uniqueidx(v) = unique(i -> v[i], eachindex(v))

# ------------------------------------------------------------
## ------------------------------------------------------------
function _histogram2D_grid(h0::Histogram, dim1, dim2;
        title = "",
        xlabel = "", ylabel = "",
        dim1_T = identity,
        dim2_T = identity,
        limits = (nothing, nothing, nothing, nothing),
        dim1_bar_width = 1.0,
        dim2_bar_width = 1.0,
        colgap = -20,
        rowgap = -5
    )

    # Plots
    f = Figure()
    g = GridLayout(f[1, 1])
    ax = Axis(g[1:3,2:5]; 
        title, limits
    )
    x1 = dim1_T(collect(keys(h0, dim1))) # koma len
    x2 = dim2_T(collect(keys(h0, dim2))) # rxn idx
    @show length(x2)
    @show cor(x1, x2)
    w = collect(values(h0))
    sidx = sortperm(w; rev = false)
    scatter!(ax, x1[sidx], x2[sidx]; 
        colormap = :viridis, markersize = 20, 
        color = log10.(w[sidx]) ./ maximum(log10, w), 
        alpha = 1.0
    )
    Colorbar(g[1:3, 6]; 
        label = "log10(count)",
        colormap = :viridis, limits = extrema(log10.(w)), 
    )
    f

    # marginals 1
    h1 = marginal(h0, dim1) # koma indx
    ax = Axis(g[4,2:5]; 
        xlabel = xlabel, ylabel = "count",
        limits = (limits[1], limits[2], nothing, nothing),
    )
    barplot!(ax, dim1_T(collect(keys(h1, 1))), collect(values(h1));
        color = :black, gap = -1, 
        width = dim1_bar_width
    )

    # marginals 2
    h1 = marginal(h0, dim2) # rxn idx
    ax = Axis(g[1:3,1]; 
        xlabel = "count", ylabel = ylabel,
        limits = (nothing, nothing, limits[3], limits[4]),
        xticklabelrotation = pi/4
    )
    barplot!(ax, dim2_T(collect(keys(h1, 1))), collect(values(h1));
        direction=:x, color = :black, gap = -1, 
        width = dim2_bar_width
    )

    colgap!(g, 1, colgap)
    rowgap!(g, 3, rowgap)
    
    f
end


# ------------------------------------------------------------
nothing