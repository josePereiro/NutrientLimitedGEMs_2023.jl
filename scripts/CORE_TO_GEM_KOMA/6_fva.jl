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
@tempcontext ["CORE_FVA" => v"0.1.0"] let
    
    # dbs
    glob_db = query(["ROOT", "GLOBALS"])
    xlep_db = query(["ROOT", "CORE_XLEP"])

    # koma files
    downreg_factor = 0.3 # TOSYNC
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
        M, N = size(core_lep0)
        
        # run
        do_save = false
        ALG_VER = context("CORE_FVA")
        info_frec = 5
        for (obji, obj) in enumerate(obj_reg)
            
            # check done
            get(obj, "core_fva.ver", :NONE) == ALG_VER && continue
            haskey(obj, "core_strip.koset") || continue
            haskey(obj, "core_feasets") || continue
            do_save = true

            # info
            show_flag = obji == 1 || obji == lastindex(obj_reg) || iszero(rem(obji, info_frec)) 
            show_flag && println("[", getpid(), ".", threadid(), "] ", 
                "obji ", obji, "\\", length(obj_reg), " ",
                basename(fn)
            )
            
            koset = obj["core_strip.koset"]
            feasets = obj["core_feasets"]
            for (li, feaobj) in feasets
                feaset = koset[1:li]

                # fva
                _with_downreg(core_lep0, feaset, downreg_factor) do
                    fvalb, fvaub = fva(core_lep0, LP_SOLVER; verbose = false)
                    feaobj["core_fva.fvalb"] = Float16.(fvalb)
                    feaobj["core_fva.fvaub"] = Float16.(fvaub)
                end 
            end
            obj["core_fva.ver"] = ALG_VER
            
        end # for reg in obj_reg
        
        # serialize
        do_save && sdat(obj_reg, fn)

    end # for fn 
end