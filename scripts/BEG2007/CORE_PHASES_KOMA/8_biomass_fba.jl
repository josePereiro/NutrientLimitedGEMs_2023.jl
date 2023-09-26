## ------------------------------------------------------------
@time begin
    using Random
    using MetXGEMs
    using MetXNetHub
    using BlobBatches
    using MetXBase
    using MetXOptim
    using NutrientLimitedGEMs
end

# ------------------------------------------------------------
include("1_setup_sim.jl")
include("1.1_utils.jl")

## ------------------------------------------------------------
@tempcontext ["CORE_BIOMASS_FBA" => v"0.1.0"] let
    
    # context
    ALG_VER = context("CORE_BIOMASS_FBA")
    GLOB_DB = query(["ROOT", "GLOBALS"])
    XLEP_DB = query(["ROOT", "CORE_XLEP"])
    LP_SOLVER = GLOB_DB["LP_SOLVER"]
    DOWNREG_FACTOR = GLOB_DB["DOWNREG_FACTOR"] 

    # read batches
    n0 = 0 # init file
    n1 = Inf # non-ignored file count
    _th_readdir(n1, n0; nthrs = 10) do bbi, bb

        # filter
        islocked(bb) && return :ignore # somebody is working
        get(bb["meta"], "core_biomass_fba.ver", :NONE) == ALG_VER && return :ignore
        haskey(bb["meta"], "core_koma.ver") || return :ignore
        haskey(bb["meta"], "core_strip.ver") || return :ignore
        haskey(bb["meta"], "core_feasets.ver") || return :ignore

        lock(bb) do
        
            # lep
            core_elep0 = XLEP_DB["core_elep0"][]
            core_lep0 = lepmodel(core_elep0)
            biom_id = extras(core_lep0, "BIOM")
            core_elep0 = nothing

            # opm
            opm = FBAOpModel(core_lep0, LP_SOLVER)

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
                koset = strip_blob["koset"]
                for (li, feaobj) in feasets_blob
                    feaset = koset[1:li]
                    # fba
                    _with_downreg(opm, feaset, DOWNREG_FACTOR) do
                        optimize!(opm)
                        feaobj["core_biomass_fba.biom"] = solution(opm, biom_id)
                        feaobj["core_biomass_fba.solution"] = Float16.(solution(opm))
                    end  # _with_downreg

                end # for feasets

                # GC
                gc_flag = iszero(rem(blobi, gc_frec))
                gc_flag && GC.gc()

            end # for blobi
            
            # sign
            bb["meta"]["core_biomass_fba.ver"] = ALG_VER
            
            # save
            serialize(bb)

        end # lock
    end # for fn 
end