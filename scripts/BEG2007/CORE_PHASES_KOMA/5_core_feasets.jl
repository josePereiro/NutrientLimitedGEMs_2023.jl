## ------------------------------------------------------------
@time begin
    using Random
    using MetXGEMs
    using MetXBase
    using MetXOptim
    using BlobBatches
    using NutrientLimitedGEMs
end

# ------------------------------------------------------------
include("1_setup_sim.jl")
include("1.1_utils.jl")

## ------------------------------------------------------------
@tempcontext ["CORE_FEASETS" => v"0.1.0"] let
    
    # context
    ALG_VER = context("CORE_FEASETS")
    GLOB_DB = query(["ROOT", "GLOBALS"])
    XLEP_DB = query(["ROOT", "CORE_XLEP"])
    DOWNREG_BATCH_SIZE = GLOB_DB["DOWNREG_BATCH_SIZE"]

    # koma files
    _th_readdir(Inf, 0; nthrs = 10) do bbi, bb

        # filter
        islocked(bb) && return :continue # somebody is working
        get(bb["meta"], "core_feasets.ver", :NONE) == ALG_VER && return :continue
        haskey(bb["meta"], "core_koma.ver") || return :continue
        haskey(bb["meta"], "core_strip.ver") || return :continue

        # lock
        lock(bb) do
            
            # new frame
            feasets_frame = bb["core_feasets"] = Dict[]

            # run
            info_frec = 100
            strip_frame = bb["core_strip"]
            for (blobi, strip_blob) in enumerate(strip_frame)
                
                # push! new obj
                feasets_blob = Dict{Int16, Any}()
                push!(feasets_frame, feasets_blob)

                # info
                info_flag = blobi == 1 || blobi == lastindex(strip_frame) || iszero(rem(blobi, info_frec)) 
                info_flag && println("[", getpid(), ".", threadid(), "] ", 
                    "blobi ", blobi, "\\", length(strip_frame), " ",
                    basename(rootdir(bb))
                )
                
                # gen feasets
                koset = strip_blob["koset"]
                idxs = range(firstindex(koset), lastindex(koset) - 1; step = DOWNREG_BATCH_SIZE)
                for lasti in idxs
                    feasets_blob[lasti] = Dict{String, Any}()
                end
            end # for reg in obj_reg
            
            # sign
            bb["meta"]["core_feasets.ver"] = ALG_VER

            # save
            serialize(bb)
            
        end # lock
    end # for bb 
end