# ------------------------------------------------------------
function _log(msg; loginfo...)

    # format log info
    ks = collect(keys(loginfo))
    # sort!(ks)
    loginfo = [string(k, "=", loginfo[k]) for k in ks]
    loginfo = string("[", getpid(), ".", threadid(), "] ", now(), " ", msg, " | ", join(loginfo, ", "))
    
    # log!
    logfn = procdir(PROJ, [SIMVER], "koma.log")
    mkpath(dirname(logfn))
    try; open((io) -> println(io, loginfo), logfn, "a"); catch ignored end
    return logfn
end

# ------------------------------------------------------------
function _sync_state!(koma_hashs, obj_reg; loginfo...)
    # save state
    _, _koma_hashs = lprocdat(PROJ, [SIMVER], "koma_hashs", ".jls") do 
        UInt64[]
    end
    push!(koma_hashs, setdiff(_koma_hashs, koma_hashs)...)
    unique!(koma_hashs)
    sort!(koma_hashs)
    
    sprocdat(PROJ, koma_hashs, 
        [SIMVER], "koma_hashs", ".jls"
    )
    sprocdat(PROJ, obj_reg, 
        [SIMVER], "obj_reg", (;h = hash(koma_hashs)), ".jls"
    )

    # log
    _log("SYNC"; h=hash(koma_hashs), loginfo...)
end


# ------------------------------------------------------------
function _foreach_obj_reg(f::Function; n = Inf)
    ch = Channel{Tuple}(nthreads(); spawn = true) do _ch
        files = readdir(procdir(PROJ, [SIMVER]); join = true)
        i = 0
        for fn in files
            contains(basename(fn), "obj_reg") || continue
            i >= n && break
            i += 1
            fn, obj_reg = ldat(fn)
            put!(_ch, (fn, obj_reg))
        end
    end
    for (fn, obj_reg) in ch
        f(fn, obj_reg) === true && break
    end
    return nothing
end

# ------------------------------------------------------------
# function _foreach_obj_reg(f::Function; n = Inf)
#     files = readdir(procdir(PROJ, [SIMVER]); join = true)
#     i = 0 
#     lk1 = ReentrantLock()
#     lk2 = ReentrantLock()
#     @threads for fn in files
#         contains(basename(fn), "obj_reg") || continue
#         lock(lk1) do
#             i >= n && return true
#             i += 1
#             return false
#         end && break
#         fn, obj_reg = ldat(fn)
#         lock(lk2) do 
#             f(fn, obj_reg) === true && return true
#             return false
#         end && break
#     end
#     return nothing
# end

function _load_and_vectorize(k::String, ks::String...; n = Inf)
    ks = String[k; ks...]
    col0 = Dict{String, Vector}(k => [] for k in ks)
    col_pool = [deepcopy(col0) for _ in 1:nthreads()]
    files = readdir(procdir(PROJ, [SIMVER]); join = true)
    lk = ReentrantLock()
    i = 0 
    @threads for fn in files
        contains(basename(fn), "obj_reg") || continue
        col = col_pool[threadid()]
        lock(lk) do
            i >= n && return true
            i += 1
            return false
        end && break
        _, obj_reg = ldat(fn)
        for k in ks  
            colv = col[k]
            for obj in obj_reg
                push!(colv, obj[k])
            end
        end
    end
    # reduce
    for col in col_pool
        for (k, vec) in col
            push!(col0[k], vec...)
        end
    end
    return col0
end

function _with_kos(f::Function, model, kos::Vector; zero = 0.0)
    kos = colindex(model, kos)
    lb0, ub0 = bounds(model, kos)
    try; bounds!(model, kos, zero, zero); f()
        finally; bounds!(model, kos, lb0, ub0)
    end
end


# ------------------------------------------------------------
nothing