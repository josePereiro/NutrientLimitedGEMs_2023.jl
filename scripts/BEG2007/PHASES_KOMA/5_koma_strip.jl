## ------------------------------------------------------------
@time begin
    using Random
    using MetXGEMs
    using MetXBase
    using MetXOptim
    using BlobBatches
    using NutrientLimitedGEMs_2023
end

# TODO: Add Graphs.jl kind of functionality for getting basic stuff
# ------------------------------------------------------------
include("1_setup.jl")
include("2_utils.jl")

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
    GLOB_DB = query(["ROOT", "GLOBALS"])
    XLEP_DB = query(["ROOT", "CORE_XLEP"])

    # koma files
    n0 = 0
    n1 = Inf
    _th_readdir(;n0, n1, nthrs = 10) do bbi, bb

        # filter
        islocked(bb) && return :continue # somebody is working
        get(bb["meta"], "core_strip.ver", :NONE) == ALG_VER && return :continue
        haskey(bb["meta"], "core_koma.ver") || return :continue

        lock(bb) do

            # new frame
            bb["core_strip"] = Dict[]

            # globals
            LP_SOLVER = GLOB_DB["LP_SOLVER"]
            
            # lep
            core_elep0 = XLEP_DB["core_elep0"][]
            core_lep0 = lepmodel(core_elep0)
            core_elep0 = nothing
            obj_id = extras(core_lep0, "BIOM")
            obj_idx = colindex(core_lep0, obj_id)

            # opm
            opm = FBAOpModel(core_lep0, LP_SOLVER)
            
            # run
            info_frec = 100
            koma_frame = bb["core_koma"]
            strip_frame = bb["core_strip"] # new frame
            for (blobi, koma_blob) in enumerate(koma_frame)
                
                # new blob
                strip_blob = typeof(koma_blob)()
                push!(strip_frame, strip_blob)

                # info
                info_flag = blobi == 1 || blobi == lastindex(koma_frame) || iszero(rem(blobi, info_frec)) 
                info_flag && println("[", getpid(), ".", threadid(), "] ", 
                    "blobi ", blobi, "\\", length(koma_frame), " ",
                    basename(rootdir(bb))
                )

                # koset
                koset = koma_blob["koset"]
                isempty(koset) && continue

                # strip
                strip_blob["koset"] = Int16[]
                DOWNREG_FACTOR = GLOB_DB["DOWNREG_FACTOR"]
                obj_val_th = GLOB_DB["KO_OBJ_VAL_TH"] 
                for downi in reverse(eachindex(koset))
                    feasible = true
                    subdown = koset[1:downi]
                    _with_downreg(opm, subdown, DOWNREG_FACTOR) do
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
        return :continue
    end # for fn 
end