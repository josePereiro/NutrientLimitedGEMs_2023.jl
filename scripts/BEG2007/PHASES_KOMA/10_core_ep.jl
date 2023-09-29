## ------------------------------------------------------------
@time begin
    using Random
    using MetXGEMs
    using MetXBase
    using MetXEP
    using BlobBatches
    using MetXOptim
    using MetXEP
    using Base.Threads
    using NutrientLimitedGEMs
end

# TODO: Add Graphs.jl kind of functionality for getting basic stuff
# ------------------------------------------------------------
include("1_setup.jl")
include("2_utils.jl")

## ------------------------------------------------------------
@tempcontext ["CORE_EP" => v"0.1.0"] let
    
    # context
    ALG_VER = context("CORE_EP")
    GLOB_DB = query(["ROOT", "GLOBALS"])
    XLEP_DB = query(["ROOT", "CORE_XLEP"])
    
    # read batches
    n0 = 0 # init file
    n1 = Inf # non-ignored file count
    _th_readdir(;n0, n1, nthrs = 10) do bbi, bb
        
        # filter
        islocked(bb) && return :ignore # somebody is working
        get(bb["meta"], "core_ep.ver", :NONE) == ALG_VER && return :ignore
        haskey(bb["meta"], "core_koma.ver") || return :ignore
        haskey(bb["meta"], "core_strip.ver") || return :ignore
        haskey(bb["meta"], "core_feasets.ver") || return :ignore
        haskey(bb["meta"], "core_fva.ver") || return :ignore

        # lock
        lock(bb) do

            # lep
            core_elep0 = XLEP_DB["core_elep0"][]
            core_lep0 = lepmodel(core_elep0)
            core_elep0 = nothing
            M, N = size(core_lep0)
            
            # run
            info_frec = 2
            gc_frec = 10
            
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
                    fvalb = get(feaobj, "core_fva.fvalb", [])
                    fvaub = get(feaobj, "core_fva.fvaub", [])
                    try
                        isempty(fvalb) && error("go to catch")
                        isempty(fvaub) && error("go to catch")
                        
                        # contextu fva bounds
                        lb!(core_lep0, fvalb)
                        ub!(core_lep0, fvaub)
                        
                        # epm
                        epm = FluxEPModelT0(core_lep0)
                        config!(epm; 
                            verbose = false,    
                            epsconv = 1e-6
                        )
                        converge!(epm)

                        feaobj["core_ep.entropy"] = entropy(epm)
                        feaobj["core_ep.free_energy"] = first(free_energy(epm))
                        feaobj["core_ep.status"] = convergence_status(epm)
                    catch e
                        feaobj["core_ep.entropy"] = NaN
                        feaobj["core_ep.free_energy"] = NaN
                        feaobj["core_ep.status"] = :error
                        # rethrow(e)
                    end
                    
                end # for (li, feaobj)
                
                # GC
                gc_flag = iszero(rem(blobi, gc_frec))
                gc_flag && GC.gc()
                
            end # for blobi
            
            # sign
            bb["meta"]["core_ep.ver"] = ALG_VER
            
            # save
            serialize(bb)
        end # lock

    end # for fn 
end