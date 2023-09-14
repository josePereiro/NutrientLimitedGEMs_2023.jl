## ------------------------------------------------------------
@time begin
    using Random
    using MetXGEMs
    using MetXBase
    using MetXOptim
    using Base.Threads
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
    glob_db = query(["ROOT", "GLOBALS"])
    xlep_db = query(["ROOT", "CORE_XLEP"])

    # koma files
    objfiles = readdir(procdir(PROJ, [SIMVER]); join = true)


    @threads for fn in shuffle(objfiles)
        contains(basename(fn), "obj_reg") || continue

        # deserialize
        _, obj_reg = ldat(fn)
        isempty(obj_reg) && continue

        # globals
        LP_SOLVER = glob_db["LP_SOLVER"]
        
        # lep
        core_elep0 = xlep_db["core_elep0"][]
        core_lep0 = lepmodel(core_elep0)
        core_elep0 = nothing
        obj_id = extras(core_lep0, "BIOM")
        obj_idx = colindex(core_lep0, obj_id)

        # run
        do_save = false
        ALG_VER = context("CORE_FEASETS")
        info_frec = 100
        for (obji, obj) in enumerate(obj_reg)

            # check done
            get(obj, "core_feasets.ver", :NONE) == ALG_VER && continue
            haskey(obj, "core_strip.koset") || continue
            do_save = true

            # info
            info_flag = obji == 1 || obji == lastindex(obj_reg) || iszero(rem(obji, info_frec)) 
            info_flag && println("[", getpid(), ".", threadid(), "] ", 
                "obji ", obji, "\\", length(obj_reg), " ",
                basename(fn)
            )
            
            # gen feasets
            koset = obj["core_strip.koset"]
            feaset_gen_step = 3 # TOSYNC
            idxs = range(firstindex(koset), lastindex(koset) - 1; step = feaset_gen_step)
            obj["core_feasets"] = Dict{Int16, Any}()
            for lasti in idxs
                obj["core_feasets"][lasti] = Dict{String, Any}()
            end
            
            obj["core_feasets.ver"] = ALG_VER

        end # for reg in obj_reg
        
        # save
        do_save && sdat(obj_reg, fn)
        
    end # for fn 
end