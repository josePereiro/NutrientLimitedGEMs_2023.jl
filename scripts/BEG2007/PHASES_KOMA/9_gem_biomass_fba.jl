## ------------------------------------------------------------
@time begin
    using Random
    using MetXGEMs
    using MetXNetHub
    using BlobBatches
    using MetXBase
    using MetXOptim
    using NutrientLimitedGEMs_2023
end

# ------------------------------------------------------------
include("1_setup.jl")
include("2_utils.jl")

## ------------------------------------------------------------
@tempcontext ["CORE_BIOMASS_FBA" => v"0.1.0"] let
    
    # context
    ALG_VER = context("CORE_BIOMASS_FBA")
    GLOB_DB = query(["ROOT", "GLOBALS"])
    CORE_XLEP_DB = query(["ROOT", "CORE_XLEP"])
    GEM_XLEP_DB = query(["ROOT", "GEM_XLEP"])
    RXN_MAP = query(["ROOT", "GEM_CORE_RXN_MAP"])["rxn_map"]
    LP_SOLVER = GLOB_DB["LP_SOLVER"]
    DOWNREG_FACTOR = GLOB_DB["DOWNREG_FACTOR"] 

    # read batches
    n0 = 0 # init file
    n1 = Inf # non-ignored file count
    _th_readdir(;n0, n1, nthrs = 10) do bbi, bb

        # filter
        islocked(bb) && return :ignore # somebody is working
        get(bb["meta"], "gem_biomass_fba.ver", :NONE) == ALG_VER && return :ignore
        haskey(bb["meta"], "core_koma.ver") || return :ignore
        haskey(bb["meta"], "core_strip.ver") || return :ignore
        haskey(bb["meta"], "core_feasets.ver") || return :ignore

        lock(bb) do
        
            # lep
            core_elep0 = CORE_XLEP_DB["core_elep0"][]
            core_lep0 = lepmodel(core_elep0)
            core_elep0 = nothing
            
            # TDOD: boxing kill GEM growth, using unbox lep0
            # gem_elep0 = GEM_XLEP_DB["gem_elep0"][]
            # gem_lep0 = lepmodel(gem_elep0)
            # gem_elep0 = nothing
            # gem_biom_id = extras(gem_lep0, "BIOM")

            gem_lep0 = GEM_XLEP_DB["gem_lep0"][]
            gem_biom_id = extras(gem_lep0, "BIOM")

            # gem_opm
            gem_opm = FBAOpModel(gem_lep0, LP_SOLVER)

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
                    core_feaset = koset[1:li]
                    gem_feaset = [RXN_MAP[core_id]
                        for core_id in colids(core_lep0, core_feaset)
                    ]

                    # fba
                    _with_downreg(gem_opm, gem_feaset, DOWNREG_FACTOR) do
                        try
                            optimize!(gem_opm)
                            feaobj["gem_biomass_fba.biom"] = solution(gem_opm, gem_biom_id)
                            feaobj["gem_biomass_fba.solution"] = Float16.(solution(gem_opm))
                        catch e
                            feaobj["gem_biomass_fba.biom"] = NaN
                            feaobj["gem_biomass_fba.solution"] = Float16[]
                        end
                    end  # _with_downreg
    
                end # for feasets

                # GC
                gc_flag = iszero(rem(blobi, gc_frec))
                gc_flag && GC.gc()

            end # for blobi
            
            # sign
            bb["meta"]["gem_biomass_fba.ver"] = ALG_VER
            
            # save
            serialize(bb)

        end # lock
    end # for fn 
end