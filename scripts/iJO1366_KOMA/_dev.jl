## ------------------------------------------------------------
@time begin
    using Dates
    using Random
    using MetXGEMs
    using MetXBase
    using MetXOptim
    using MetXNetHub
    using Base.Threads
    using ProgressMeter
    using SimpleLockFiles
    using NutrientLimitedGEMs
end

# ------------------------------------------------------------
include("1_setup_sim.jl")
include("1.1_utils.jl")
 
## ------------------------------------------------------------
# TODO: sync ECOLI_CORE and iJO1366 code
# DONE: start running the trimming script
# DONE: Finish the review
# DONE: Move to cluster SingleCell_Lactate_2023
let
    files = readdir(procdir(PROJ, [SIMVER]); join = true)
    for fn in files
        startswith(basename(fn), "obj_reg") || continue
        _, obj_reg = ldat(fn)
        do_save = false
        for obj in obj_reg
            haskey(obj, "strip.koset") || continue
            if isempty(obj["strip.koset"])
                obj["strip.koset"] = Int16[]
                do_save = true
            end
        end
        do_save || continue
        println("[", getpid(), ".", threadid(), "] ", fn)
        sdat(obj_reg, fn)
    end
end