## ------------------------------------------------------------
@time begin
    using Random
    using MetXGEMs
    using MetXBase
    using MetXOptim
    using BlobBatches
    using NutrientLimitedGEMs
end

# TODO: Add Graphs.jl kind of functionality for getting basic stuff
# ------------------------------------------------------------
include("1_setup_sim.jl")
include("1.1_utils.jl")

## ------------------------------------------------------------
# IDEA: ProjFlows, save all variables except the one with a given naming convention 
# Ex: var (yes), _var (yes), __var (no)
# Also, have a gitignore kind of configuration
# All computation/saving happends at a @commit like macro

# TODO: Add file lock to 'sdat'...

# IDEA: Add obj renamer cli

@tempcontext ["CORE_FEASETS" => v"0.1.0"] let
    
    # dbs
    ALG_VER = context("CORE_FEASETS")
    glob_db = query(["ROOT", "GLOBALS"])
    xlep_db = query(["ROOT", "CORE_XLEP"])

    # koma files
    batches = readdir(BlobBatch, procdir(PROJ, [SIMVER]))
    for (bbi, bb) in enumerate(batches)

        # filter
        islocked(bb) && continue # somebody is working
        get(bb["meta"], "core_feasets.ver", :NONE) == ALG_VER && continue
        haskey(bb["meta"], "core_koma.ver") || continue
        haskey(bb["meta"], "core_strip.ver") || continue

        # lock
        lock(bb) do
            
            # new frame
            bb["core_feasets"] = Dict[]

            # run
            info_frec = 100
            for (blobi, strip_blob) in enumerate(bb["core_strip"])
                
                # feasibles
                feasets_blob = typeof(strip_blob)()
                push!(bb["core_feasets"], feasets_blob)

                # info
                info_flag = blobi == 1 || blobi == lastindex(bb["core_strip"]) || iszero(rem(blobi, info_frec)) 
                info_flag && println("[", getpid(), ".", threadid(), "] ", 
                    "blobi ", blobi, "\\", length(bb["core_strip"]), " ",
                    basename(rootdir(bb))
                )
                
                # gen feasets
                koset = strip_blob["koset"]
                feaset_gen_step = 3 # TOSYNC
                idxs = range(firstindex(koset), lastindex(koset) - 1; step = feaset_gen_step)
                feasets_blob["feasets"] = Dict{Int16, Any}()
                for lasti in idxs
                    feasets_blob["feasets"][lasti] = Dict{String, Any}()
                end

            end # for reg in obj_reg
            
            # sign
            bb["meta"]["core_feasets.ver"] = ALG_VER

            # save
            serialize(bb)
            
        end # lock
    end # for bb 
end