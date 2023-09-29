## --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# Collect max biomass fba solutions
let
    # context
    _simver = "ECOLI-CORE-BEG2007-PHASE_I-0.1.0"
    _load_contextdb(_simver)

    n0 = 0 # init file
    n1 = Inf # non-ignored file count
    cid = (@__FILE__, _simver, "ep:entropy + free energy", n0, n1)
    lk = ReentrantLock()
    _, ret = withcachedat(PROJ, :get!, cid) do
        _h0 = Histogram(
                   0:10000,                        # _fealen
            -10000.0:0.1:10000.0,                  # entropy
            -10000.0:0.1:10000.0,                  # free energy
        )
        h_pool = Dict()
        _th_readdir(_simver; n0, n1, nthrs = 10) do bbi, bb
            haskey(bb["meta"], "core_ep.ver") || return :ignore
            # load frame
            feasets_db = bb["core_feasets"]
            _h = get!(() -> deepcopy(_h0), h_pool, threadid())
            for feasets_blob0 in feasets_db
                for (_fealen, feasets_blob1) in feasets_blob0
                    feaobj["core_ep.status"] == :error && continue
                    _H = feaobj["core_ep.entropy"]
                    _F = feaobj["core_ep.free_energy"]
                    count!(_h, (_fealen, _H, _F))
                end
            end # for feasets_blob0
        end # _th_readdir
        # reduce
        merge!(_h0, h_pool...) # reduce
        return _h0
    end
    global h) = ret

    return
end