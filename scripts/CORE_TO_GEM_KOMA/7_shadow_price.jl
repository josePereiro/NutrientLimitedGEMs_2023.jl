## ------------------------------------------------------------
@time begin
    using Random
    using MetXGEMs
    using MetXNetHub
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
@tempcontext ["CORE_NUT_SP" => v"0.1.0"] let
    
    # dbs
    ALG_VER = context("CORE_NUT_SP")
    glob_db = query(["ROOT", "GLOBALS"])
    xlep_db = query(["ROOT", "CORE_XLEP"])

    # koma files
    downreg_factor = 0.3 # TOSYNC
    objfiles = readdir(procdir(PROJ, [SIMVER]); join = true)
    for fn in shuffle(objfiles)
        contains(basename(fn), "obj_reg") || continue

        # deserialize
        obj_reg = try_ldat(fn)
        isempty(obj_reg) && continue
        
        # globals
        LP_SOLVER = glob_db["LP_SOLVER"]
        
        # lep
        core_elep0 = xlep_db["core_elep0"][]
        core_lep0 = lepmodel(core_elep0)
        core_elep0 = nothing

        # opm
        opm = FBAOpModel(core_lep0, LP_SOLVER)
        
        # run
        do_save = false
        info_frec = 2
        gc_frec = 10
        for (obji, obj) in enumerate(obj_reg)
            
            # check done
            get(obj, "core_nut_sp.ver", :NONE) == ALG_VER && continue
            haskey(obj, "core_strip.koset") || continue
            haskey(obj, "core_feasets") || continue
            do_save = true
            
            # info
            show_flag = obji == 1 || obji == lastindex(obj_reg) || iszero(rem(obji, info_frec)) 
            show_flag && println("[", getpid(), ".", threadid(), "] ", 
                "obji ", obji, "\\", length(obj_reg), " ",
                basename(fn)
            )
            
            # compute
            koset = obj["core_strip.koset"]
            feasets = obj["core_feasets"]
            for (li, feaobj) in feasets
                feaset = koset[1:li]
                _with_downreg(opm, feaset, downreg_factor) do

                    for exch in [
                            "EX_glc__D_e", "EX_lac__D_e", "EX_malt_e",
                            "EX_gal_e", "EX_glyc_e", "EX_ac_e", "EX_ac_e"
                        ]

                        test_points = lb(opm, exch) .* [1.0, 0.95, 0.9]
                        obj_m, obj_err, vars_ms, vars_errs = bound_dual_prices(
                            opm, exch, test_points, :lb; 
                            dovars = true
                        )
                        feaobj["core_nut_sp.$(exch).obj_m"] = Float16.(obj_m)
                        feaobj["core_nut_sp.$(exch).obj_err"] = Float16.(obj_err)
                        feaobj["core_nut_sp.$(exch).vars_ms"] = Float16.(vars_ms)
                        feaobj["core_nut_sp.$(exch).vars_errs"] = Float16.(vars_errs)
                        
                    end # for exch
                end  # _with_downreg
            end # for feasets
            
            # ALG_VER
            obj["core_nut_sp.ver"] = ALG_VER

            # GC (TEST)
            gc_flag = iszero(rem(obji, gc_frec))
            gc_flag && GC.gc()
            
        end # for reg in obj_reg
        
        # serialize
        do_save && sdat(obj_reg, fn)

        # TEST
        do_save && exit()

    end # for fn 
end