## ------------------------------------------------------------
@time begin
    using MetXGEMs
    using MetXNetHub
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
# TODD: Update to BlobBatches
@tempcontext ["GEM_NUT_SP" => v"0.1.0"] let
    
    # context
    ALG_VER = context("CORE_NUT_SP")
    GLOB_DB = query(["ROOT", "GLOBALS"])
    CORE_XLEP_DB = query(["ROOT", "CORE_XLEP"])
    GEM_XLEP_DB = query(["ROOT", "GEM_XLEP"])
    RXN_MAP = query(["ROOT", "GEM_CORE_RXN_MAP"])["rxn_map"]
    LP_SOLVER = GLOB_DB["LP_SOLVER"]
    DOWNREG_FACTOR = GLOB_DB["DOWNREG_FACTOR"] # TOSYNC

    # read batches
    n0 = 0 # init file
    n1 = Inf # non-ignored file count
    _th_readdir(;n0, n1, nthrs = 10) do bbi, bb

        # filter
        islocked(bb) && return :ignore # somebody is working
        get(bb["meta"], "gem_nut_sp.ver", :NONE) == ALG_VER && return :ignore
        haskey(bb["meta"], "core_koma.ver") || return :ignore
        haskey(bb["meta"], "core_strip.ver") || return :ignore
        haskey(bb["meta"], "core_feasets.ver") || return :ignore
        
        lock(bb) do

            # lep
            core_elep0 = CORE_XLEP_DB["core_elep0"][]
            core_lep0 = lepmodel(core_elep0)
            core_elep0 = nothing

            gem_elep0 = GEM_XLEP_DB["core_elep0"][]
            gem_lep0 = lepmodel(gem_elep0)
            gem_elep0 = nothing

            # gem_opm
            gem_opm = FBAOpModel(gem_elep0, LP_SOLVER)

            # frames
            feasets_frame = bb["core_feasets"]
            strip_frame = bb["core_strip"]

            # run
            info_frec = 50
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
                core_koset = strip_blob["koset"]
                for (li, feaobj) in feasets_blob
                    core_feaset = core_koset[1:li]
                    gem_feaset = [RXN_MAP[core_id]
                        for core_id in colids(core_lep0, core_feaset)
                    ]

                    _with_downreg(gem_opm, gem_feaset, DOWNREG_FACTOR) do
                        for exch in [
                                "EX_glc__D_e", "EX_lac__D_e", "EX_malt_e",
                                "EX_gal_e", "EX_glyc_e", "EX_ac_e"
                            ]
                            exch in colids(gem_opm) || continue
                            test_points = lb(gem_opm, exch) .* [1.0, 0.95, 0.9]
                            obj_m, obj_err, vars_ms, vars_errs = try
                                bound_dual_prices(
                                    gem_opm, exch, test_points, :lb; 
                                    dovars = true
                                )
                            catch e
                                NaN, NaN, [], []
                            end
                            feaobj["gem_nut_sp.$(exch).obj_m"] = Float16(obj_m)
                            feaobj["gem_nut_sp.$(exch).obj_err"] = Float16(obj_err)
                            feaobj["gem_nut_sp.$(exch).vars_ms"] = Float16.(vars_ms)
                            feaobj["gem_nut_sp.$(exch).vars_errs"] = Float16.(vars_errs)
                            
                        end # for exch
                    end  # _with_downreg

                end # for feasets
                
                # GC
                gc_flag = iszero(rem(blobi, gc_frec))
                gc_flag && GC.gc()
                
            end # for blobi

            # sign
            bb["meta"]["gem_nut_sp.ver"] = ALG_VER
            
            # save
            serialize(bb)
        end

    end # for fn 
end