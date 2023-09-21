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

@tempcontext ["CORE_STRIP" => v"0.1.0"] let
    
    # dbs
    ALG_VER = context("CORE_STRIP")
    glob_db = query(["ROOT", "GLOBALS"])
    xlep_db = query(["ROOT", "CORE_XLEP"])

    # koma files
    batches = readdir(BlobBatch, procdir(PROJ, [SIMVER]))
    for (bbi, bb) in enumerate(batches)

        # filter
        islocked(bb) && continue # somebody is working
        get(bb["meta"], "core_strip.ver", :NONE) == ALG_VER && continue
        haskey(bb["meta"], "core_koma.ver") || continue

        # lock
        lock(bb) do

            # new frame
            bb["core_strip"] = Dict[]

            # globals
            LP_SOLVER = glob_db["LP_SOLVER"]
            
            # lep
            core_elep0 = xlep_db["core_elep0"][]
            core_lep0 = lepmodel(core_elep0)
            core_elep0 = nothing
            obj_id = extras(core_lep0, "BIOM")
            obj_idx = colindex(core_lep0, obj_id)

            # opm
            opm = FBAOpModel(core_lep0, LP_SOLVER)
            
            # run
            info_frec = 100
            for (blobi, koma_blob) in enumerate(bb["core_koma"])
                
                # push blob
                strip_blob = typeof(koma_blob)()
                push!(bb["core_strip"], strip_blob)

                # info
                info_flag = blobi == 1 || blobi == lastindex(bb["core_koma"]) || iszero(rem(blobi, info_frec)) 
                info_flag && println("[", getpid(), ".", threadid(), "] ", 
                    "blobi ", blobi, "\\", length(bb["core_koma"]), " ",
                    basename(rootdir(bb))
                )

                # koset
                koset = koma_blob["koset"]
                isempty(koset) && continue

                # strip
                strip_blob["koset"] = Int16[]
                downreg_factor = 0.3 # TOSYNC
                obj_val_th = 0.01 # TOSYNC
                for downi in reverse(eachindex(koset))
                    feasible = true
                    subdown = koset[1:downi]
                    _with_downreg(opm, subdown, downreg_factor) do
                        set_linear_obj!(opm, obj_idx, MAX_SENSE)
                        sol = 0.0
                        try; optimize!(opm)
                            sol = solution(opm, obj_idx)
                        catch e; end
                        feasible = sol > obj_val_th
                    end
                    
                    if feasible 
                        strip_blob["koset"] = koset[1:downi+1]
                        break
                    end 
                end # for downi 
            end # for core_blob

            # sign
            bb["meta"]["core_strip.ver"] = ALG_VER
            
            # save
            serialize(bb)
            
        end # lock
    end # for fn 
end