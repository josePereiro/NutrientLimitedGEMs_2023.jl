## ------------------------------------------------------------
import Base.Threads: Atomic
function _th_readdir(f::Function, sim; 
        n1 = Inf, n0 = 1,
        info_frec = 1.0, nthrs = 10, verbose = true, 
        perm = sort!
    )
    batches = readdir(BlobBatch, procdir(PROJ, [sim]); perm)
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
import Optim
function _MaxEnt_beta(av0)
    a0, b0 = 1.0, 1.0
    a, b = 1.0, 1.0
    res_b = Optim.optimize(0.0, Inf, [b0]) do _b
        b = clamp(first(_b), 1e-9, Inf)
        res_a = Optim.optimize(0.0, Inf, [a0]) do _a
            a = clamp(first(_a), 1e-9, Inf)
            B = Beta(a, b)
            (av0 - mean(B))^2
        end
        _a = first(Optim.minimizer(res_a))
        a = clamp(first(_a), 1e-9, Inf)
        B = Beta(a, b)
        S = entropy(B)
        return -S
    end
    b = first(Optim.minimizer(res_b))
    B = Beta(a, b)
    return B
end

# ------------------------------------------------------------
nothing