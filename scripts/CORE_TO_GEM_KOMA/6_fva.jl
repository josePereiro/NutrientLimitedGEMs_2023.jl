## ------------------------------------------------------------
@time begin
    using Random
    using MetXGEMs
    using MetXBase
    using BlobBatches
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
    ALG_VER = context("CORE_FVA")
    glob_db = query(["ROOT", "GLOBALS"])
    xlep_db = query(["ROOT", "CORE_XLEP"])
    LP_SOLVER = glob_db["LP_SOLVER"]

    # downregulation
    downreg_factor = 0.3 # TOSYNC
    
    # read batches
    n0 = 0 # init file
    n1 = Inf # non-ignored file count
    # TODO: _th_readdir(n1, n0) makes no sense
    _th_readdir(n1, n0; nthrs = 10) do bbi, bb
        
        # filter
        islocked(bb) && return :ignore # somebody is working
        if get(bb["meta"], "core_fva.ver", :NONE) == ALG_VER 
            println("[", getpid(), ".", threadid(), "] IGNORED ", "bbi ", bbi, " ",
                "core_fva.ver ", ALG_VER, " ",
                basename(rootdir(bb))
            )
            return :ignore
        end
        haskey(bb["meta"], "core_koma.ver") || return :ignore
        haskey(bb["meta"], "core_strip.ver") || return :ignore
        haskey(bb["meta"], "core_feasets.ver") || return :ignore
        haskey(bb["meta"], "core_biomass.ver") || return :ignore
        haskey(bb["meta"], "core_nut_sp.ver") || return :ignore

        # lock
        # lock(bb) do

            # lep
            core_elep0 = xlep_db["core_elep0"][]
            core_lep0 = lepmodel(core_elep0)
            core_elep0 = nothing
            M, N = size(core_lep0)
            
            # run
            info_frec = 2
            
            # frames
            feasets_frame = bb["core_feasets"]
            strip_frame = bb["core_strip"]
            
            for (blobi, (strip_blob, feasets_blob)) in enumerate(
                    zip(strip_frame, feasets_frame)
                )

                # info
                info_flag = blobi == 1 || blobi == lastindex(feasets_frame) || iszero(rem(blobi, info_frec)) 
                info_flag && println("[", getpid(), ".", threadid(), "] ", 
                    "blobi ", blobi, "\\", length(feasets_frame), " ",
                    basename(rootdir(bb))
                )
                
                koset = strip_blob["koset"]
                for (li, feaobj) in feasets_blob
                    feaset = koset[1:li]

                    # fva
                    _with_downreg(core_lep0, feaset, downreg_factor) do
                        fvalb, fvaub = fva(core_lep0, LP_SOLVER; verbose = false)
                        feaobj["core_fva.fvalb"] = Float16.(fvalb)
                        feaobj["core_fva.fvaub"] = Float16.(fvaub)
                    end 
                end
                
            end # for obj in frame
            
            # sign
            bb["meta"]["core_fva.ver"] = ALG_VER
            
            # save
            serialize(bb)
        # end # lock

    end # for fn 
end