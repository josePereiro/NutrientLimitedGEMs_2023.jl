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
using Base.Threads: Atomic
function _th_readdir(f::Function, n1 = Inf, n0 = 1;
        info_frec = 1.0, nthrs = 10, verbose = true
    )
    batches = readdir(BlobBatch, procdir(PROJ, [SIMVER]))
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

# ------------------------------------------------------------
_uniqueidx(v) = unique(i -> v[i], eachindex(v))

# ------------------------------------------------------------
nothing