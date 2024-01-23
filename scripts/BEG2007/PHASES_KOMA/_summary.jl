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
include("1_setup.jl")
include("2_utils.jl")

## ------------------------------------------------------------
# BlobBatch lock as an kwargs... serialize/deserialize
let
    n0 = 0 # init
    n1 = Inf # non-ignored count
    lk = ReentrantLock()
    while true

        try
            meta_records = Dict()
            for _SIMVER in [
                    "ECOLI-CORE-BEG2007-PHASE_0",
                    "ECOLI-CORE-BEG2007-PHASE_1",
                    "ECOLI-CORE-BEG2007-PHASE_2",
                    "ECOLI-CORE-BEG2007-PHASE_3",
                ]
                record = get!(meta_records, _SIMVER) do
                    Dict()
                end
                _th_readdir(_SIMVER; n1, n0, nthrs = 10, verbose = false) do bbi, bb
                    lock(lk) do
                        for key in keys(bb["meta"])
                            get!(record, key, 0)
                            record[key] += 1
                        end     
                        get!(record, "batches", 0)
                        record["batches"] += 1
                    end
                end
            end # for _SIMVER
            
            # Printl
            run(`clear`)
            println("="^60)
            for (ver, records) in meta_records
                println("-"^60)
                println("SIMVER: ", ver)
                batches = records["batches"]
                for (k, c) in records
                    println(k, ": ", c, " -- ", round(c / batches; sigdigits = 2))
                end
            end

            sleep(60)
        catch e
            e isa InterruptException && exit()
            rethrow(e)
        end
    end # while true
end

