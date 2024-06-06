## ------------------------------------------------------------
@time begin
    using Random
    using MetXGEMs
    using MetXBase
    using MetXOptim
    using Statistics
    using MetXNetHub
    using CairoMakie
    using Base.Threads
    using ProgressMeter
    using Combinatorics
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

@tempcontext ["STRIP" => v"0.1.0"] let
    
    # dbs
    glob_db = query(["ROOT", "GLOBALS"])
    xlep_db = query(["ROOT", "XLEP"])

    # koma files
    objfiles = readdir(procdir(PROJ, [SIMVER]); join = true)
    @threads :static for fn in shuffle(objfiles)
        contains(basename(fn), "obj_reg") || continue
        println("[", getpid(), ".", threadid(), "] ", fn)

        # deserialize
        _, obj_reg = ldat(fn)

        # globals
        LP_SOLVER = glob_db["LP_SOLVER"]
        
        # lep
        lep0 = lepmodel(xlep_db["elep0"][])
        obj_id = extras(lep0, "BIOM")
        obj_idx = colindex(lep0, obj_id)
        
        # run
        do_save = false
        alg_ver = context("STRIP")
        for obj in obj_reg

            # check done
            get(obj, "strip_ver", :NONE) == alg_ver && break
            do_save = true
            
            koset = obj["koset"]
            
            if isempty(koset) 
                obj["strip_ver"] = alg_ver
                continue
            end

            # strip
            obj["strip.koset"] = Int16[]
            for koi in reverse(eachindex(koset))
                
                ko_th = 1e-5
                feasible = true
                _with_kos(lep0, koset[1:koi]) do
                    opm = FBAOpModel(lep0, LP_SOLVER)
                    set_linear_obj!(opm, obj_idx, MAX_SENSE)
                    sol = 0.0
                    try; optimize!(opm)
                        sol = solution(opm, obj_idx)
                    catch e; end
                    feasible = sol > ko_th
                end
                
                if feasible 
                    obj["strip.koset"] = koset[1:koi+1]
                    break
                end
                
            end
            obj["strip_ver"] = alg_ver
        end # for reg in obj_reg
        
        # serialize
        do_save && sdat(obj_reg, fn)
        
    end # for fn 
end