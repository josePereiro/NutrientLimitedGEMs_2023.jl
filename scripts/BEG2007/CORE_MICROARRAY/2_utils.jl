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

# ------------------------------------------------------------
_uniqueidx(v) = unique(i -> v[i], eachindex(v))

# ------------------------------------------------------------
nothing