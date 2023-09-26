## ------------------------------------------------------------
@time begin
    using Random
    using MetXGEMs
    using MetXNetHub
    using MetXBase
    using MetXOptim
    using BlobBatches
    using BlobBatches: loadallframe!
    using Base.Threads
    using NutrientLimitedGEMs
end

# ------------------------------------------------------------
include("1_setup_sim.jl")
include("1.1_utils.jl")

## ------------------------------------------------------------
# BlobBatch lock as an kwargs... serialize/deserialize
let
    n0 = 0 # init
    n1 = Inf # non-ignored count
    lk = ReentrantLock()
    while true
        meta_count = Dict()
        for _SIMVER in [
                "ECOLI-CORE-BEG2007-PHASE_I-0.1.0",
                "ECOLI-CORE-BEG2007-PHASE_II-0.1.0",
                "ECOLI-CORE-BEG2007-PHASE_III-0.1.0",
            ]
            global SIMVER = _SIMVER
            simdat = get!(meta_count, _SIMVER) do
                Dict()
            end
            _th_readdir(n1, n0; nthrs = 10) do bbi, bb
                lock(lk) do
                    for key in keys(bb["meta"])
                        get!(simdat, key, 0)
                        simdat[key] += 1
                    end     
                    get!(simdat, "batches", 0)
                    simdat["batches"] += 1
                end
            end
        end
        
        # Printl
        for (ver, keys) in meta_count
            println("-"^60)
            println("SIMVER: ", ver)
            for (k, c) in keys
                println(k, ": ", c)
            end
        end

        sleep(60)
    end # while true
end

