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
# TODD: Update to BlobBatches
@tempcontext ["CORE_NUT_SP" => v"0.1.0"] let
    
    # dbs
    ALG_VER = context("CORE_NUT_SP")
    glob_db = query(["ROOT", "GLOBALS"])
    xlep_db = query(["ROOT", "CORE_XLEP"])
    LP_SOLVER = glob_db["LP_SOLVER"]

    # read batches
    n0 = 0 # init file
    n1 = Inf # non-ignored file count
    _th_readdir(n1, n0; nthrs = 10) do bbi, bb

        # filter
        islocked(bb) && return :ignore # somebody is working
        get(bb["meta"], "core_nut_sp.ver", :NONE) == ALG_VER && return :ignore
        haskey(bb["meta"], "core_koma.ver") || return :ignore
        haskey(bb["meta"], "core_strip.ver") || return :ignore
        haskey(bb["meta"], "core_feasets.ver") || return :ignore
        
        lock(bb) do

            # lep
            core_elep0 = xlep_db["core_elep0"][]
            core_lep0 = lepmodel(core_elep0)
            core_elep0 = nothing

            # opm
            opm = FBAOpModel(core_lep0, LP_SOLVER)

            # frames
            feasets_frame = bb["core_feasets"]
            strip_frame = bb["core_strip"]

            # run
            info_frec = 100
            gc_frec = 100
            for (blobi, (strip_blob, feasets_blob)) in enumerate(
                    zip(strip_frame, feasets_frame)
                )
                
                # info
                info_flag = blobi == 1 || blobi == lastindex(feasets_frame) || iszero(rem(blobi, info_frec)) 
                info_flag && println("[", getpid(), ".", threadid(), "] ", 
                    "blobi ", blobi, "\\", length(feasets_frame), " ",
                    basename(rootdir(bb))
                )
                
                # compute
                koset = strip_blob["koset"]
                for (li, feaobj) in feasets_blob
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
                            feaobj["core_nut_sp.$(exch).obj_m"] = Float16(obj_m)
                            feaobj["core_nut_sp.$(exch).obj_err"] = Float16(obj_err)
                            feaobj["core_nut_sp.$(exch).vars_ms"] = Float16.(vars_ms)
                            feaobj["core_nut_sp.$(exch).vars_errs"] = Float16.(vars_errs)
                            
                        end # for exch
                    end  # _with_downreg
                end # for feasets
                
                # GC (TEST)
                gc_flag = iszero(rem(obji, gc_frec))
                gc_flag && GC.gc()
                
            end # for blobi

            # sign
            bb["meta"]["core_nut_sp.ver"] = ALG_VER
            
            # save
            serialize(bb)
        end

    end # for fn 
end