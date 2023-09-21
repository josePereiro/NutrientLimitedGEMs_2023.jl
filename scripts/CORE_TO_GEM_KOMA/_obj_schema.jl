# ## ------------------------------------------------------------
# @time begin
#     using Random
#     using MetXGEMs
#     using MetXNetHub
#     using MetXBase
#     using MetXOptim
#     using BlobBatches
#     using BlobBatches: loadallframe!
#     using Base.Threads
#     using NutrientLimitedGEMs
# end

# # TODO: Add Graphs.jl kind of functionality for getting basic stuff
# # --------------------------------------------------------------
# include("1_setup_sim.jl")
# include("1.1_utils.jl")

# ## ------------------------------------------------------------- 
# # kill 1686047 1649702 1649765 1649824 1650092 1684824 1685066 1649973 1662998 1649884 1650747 1650713
# ## -------------------------------------------------------------
# function _create_bbs(fn)
#     println("="^50)
#     @show fn
#     @time obj_reg = try_ldat(fn)
#     isempty(obj_reg) && return
    
#     bb = nothing
#     for (obji, obj) in enumerate(obj_reg)
#         iszero(rem(obji, 100)) && println("obji: ", obji, " ", basename(fn))
        
#         if isnothing(bb)
#             bb = BlobBatch(procdir(PROJ, [SIMVER], "batch", (;h=hash((fn, obji)))))
#             rm(bb; force = true)
#             bb["core_koma"] = Dict[]
#             bb["core_strip"] = Dict[]
#             bb["core_feasets"] = Dict[]
#         end

#         # core_koma
#         if haskey(obj, "core_koma.ver")
#             # "core_koma.ver"    => v"0.1.0"
#             # "core_koma.status" => :ERROR
#             # "core_koma.koset"  => Int16[...]
#             bb["meta"]["core_koma.ver"] = obj["core_koma.ver"]
#             blob = Dict(
#                 "koset" => Int16.(obj["core_koma.koset"]),
#                 "status" => obj["core_koma.status"]
#             )
#             push!(bb["core_koma"], blob)
#         end

#         # core_strip
#         if haskey(obj, "core_strip.ver")
#             # "core_strip.ver"    => v"0.1.0"
#             # "core_strip.koset"  => Int16[...]
#             bb["meta"]["core_strip.ver"] = obj["core_strip.ver"]
#             blob = Dict(
#                 "koset" => Int16.(obj["core_strip.koset"]),
#             )
#             push!(bb["core_strip"], blob)
#         end

#         if haskey(obj, "core_feasets.ver")
#             bb["meta"]["core_feasets.ver"] = obj["core_feasets.ver"]
#             bb["meta"]["core_nut_sp.ver"] = obj["core_nut_sp.ver"]
#             bb["meta"]["core_biomass.ver"] = obj["core_biomass.ver"]
#             blob = obj["core_feasets"]
#             for (k, dat) in blob
#                 dat isa Vector || continue
#                 blob[k] = Float16.(dat)
#             end
#             push!(bb["core_feasets"], blob)
#         end

#         # split
#         if length(bb["core_koma"]) == 1000
#             println("-"^40)
#             @time serialize(bb)
#             @show obji
#             println(bb)
#             bb = nothing # reset
#         end
#     end
#     return nothing
# end

# ## -------------------------------------------------------------
# # I'll separate the batch files into single obje files
# let
#     bashfiles = readdir(procdir(PROJ, [SIMVER]); join = true)
#     @threads for fn in sort(bashfiles)
#         contains(basename(fn), "obj_reg") || continue
#         try
#             _create_bbs(fn)
#             rm(fn; force = true)
#         catch e
#             @error(fn)
#         end
#     end
# end

# ## -------------------------------------------------------------
# # let
# #     batches = readdir(BlobBatch, procdir(PROJ, [SIMVER]))
# #     for (bbi, bb) in enumerate(batches)
# #         global _bb = bb
# #         # @show isempty(bb)
# #         # # @show rm(bb)
# #         # @show mkpath(bb)
# #         # @show isempty(bb)
# #         # @show rm(bb)
# #     end
# # end

# ## -------------------------------------------------------------
# # let
# #     lenv = Float64[]
# #     obj_mv = Float64[]
# #     exchid = "EX_ac_e" 
# #     # exchid = "EX_lac__D_e" 
# #     # exchid = "EX_glyc_e" 
# #     # exchid = "EX_gal_e" 
# #     # exchid = "EX_glc__D_e" 
# #     # exchid = "EX_malt_e" 
# #     for blob0 in _bb["core_feasets"]
# #         for (li, blob1) in blob0
# #             push!(lenv, li)
# #             m_key = string("core_nut_sp.", exchid, ".obj_m")
# #             push!(obj_mv, blob1[m_key])
# #         end
# #     end
# #     f = Figure()
# #     ax = Axis(f[1,1]; title = exchid, xlabel = "downreg steps", ylabel = "shadown price")
# #     scatter!(ax, lenv, obj_mv)
# #     f
# # end

# ## -------------------------------------------------------------
# # ## -------------------------------------------------------------
# # let
# #     bb = BlobBatch("/home/pereiro/dev/NutrientLimitedGEMs/data/processed/CORE_TO_GEM-ECOLI-0.1.0/batch...<<h=9614697493927991082>>")
# #     loadallframe!(bb)
# #     bb
# # end
# # ## -------------------------------------------------------------
# # ## -------------------------------------------------------------
# # ## -------------------------------------------------------------
# # let
# #     # koma files
# #     downreg_factor = 0.3 # TOSYNC
# #     objfiles = readdir(procdir(PROJ, [SIMVER]); join = true)
# #     # @threads 
# #     for fn in sort(objfiles)
# #         contains(basename(fn), "obj_reg") || continue

# #         # deserialize
# #         obj_reg = try_ldat(fn)
# #         isempty(obj_reg) && continue
        
# #         info_frec = 1000
# #         do_save = false
# #         for (obji, obj) in enumerate(obj_reg)

# #             # filter
# #             haskey(obj, "core_fva.ver") || continue
# #             haskey(obj, "core_biomass.ver") || continue
# #             haskey(obj, "core_nut_sp.ver") || continue
# #             haskey(obj, "core_feasets.ver") || continue
# #             haskey(obj, "core_feasets") || continue

# #             return obj

# #             # info
# #             show_flag = obji == 1 || obji == lastindex(obj_reg) || iszero(rem(obji, info_frec)) 
# #             show_flag && println("[", getpid(), ".", threadid(), "] ", 
# #                 "obji ", obji, "\\", length(obj_reg), " ",
# #                 basename(fn)
# #             )

# #             # format
# #             # DEPRECATE

# #             # DELETE
# #             for (li, feaobj) in obj["feasibles"]
# #                 for key in feaobj
# #                     endswith(key, "vars_ms") || endswith(key, "vars_errs") || continue
# #                     feaobj[key] isa Vector{Float16} && continue
# #                     feaobj[key] = Float16.(feaobj[key])
# #                     do_save = true
# #                 end
# #             end

# #         end # for obj
        
# #         # return obj_reg
# #         do_save && sdat(obj_reg, fn)

# #     end # for fn 

# # end