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

@tempcontext ["CORE_STRIP" => v"0.1.0"] let
    
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

        # opm
        opm = FBAOpModel(core_lep0, LP_SOLVER)
        
        # run
        do_save = false
        ALG_VER = context("CORE_STRIP")
        info_frec = 100
        for (obji, obj) in enumerate(obj_reg)

            # check done
            get(obj, "core_strip.ver", :NONE) == ALG_VER && continue
            do_save = true

            # info
            info_flag = obji == 1 || obji == lastindex(obj_reg) || iszero(rem(obji, info_frec)) 
            info_flag && println("[", getpid(), ".", threadid(), "] ", 
                "obji ", obji, "\\", length(obj_reg), " ",
                basename(fn)
            )
            
            # koset
            koset = obj["core_koma.koset"]
            if isempty(koset) 
                obj["core_strip.ver"] = ALG_VER
                continue
            end

            # strip
            obj["core_strip.koset"] = Int16[]
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
                    obj["core_strip.koset"] = koset[1:downi+1]
                    break
                end 
            end # for downi 
            obj["core_strip.ver"] = ALG_VER
        end # for reg in obj_reg
        
        # save
        do_save && sdat(obj_reg, fn)
        
    end # for fn 
end