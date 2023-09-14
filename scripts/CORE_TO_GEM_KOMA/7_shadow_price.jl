## ------------------------------------------------------------
@time begin
    using Random
    using MetXGEMs
    using MetXNetHub
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
@tempcontext ["CORE_NUT_SP" => v"0.1.0"] let
    
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

        # opm
        opm = FBAOpModel(core_lep0, LP_SOLVER)
        
        # run
        do_save = false
        ALG_VER = context("CORE_NUT_SP")
        info_frec = 10
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
            
            koset = obj["core_strip.koset"]
            feaset = koset[1:(end-1)]
            
            if isempty(feaset) 
                obj["core_nut_sp.ver"] = ALG_VER
                continue
            end

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
                        feaobj["core_nut_sp.$(exch).obj_m"] = obj_m
                        feaobj["core_nut_sp.$(exch).obj_err"] = obj_err
                        feaobj["core_nut_sp.$(exch).vars_ms"] = vars_ms
                        feaobj["core_nut_sp.$(exch).vars_errs"] = vars_errs
                        
                    end # for exch
                end  # _with_downreg
            end # for feasets

            obj["core_nut_sp.ver"] = ALG_VER
            
        end # for reg in obj_reg
        
        # serialize
        do_save && sdat(obj_reg, fn)

    end # for fn 
end