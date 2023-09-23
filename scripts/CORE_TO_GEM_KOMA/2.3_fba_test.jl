# # -------------------------------------------
# # FBA TEST

# println()
# println("="^40)
# println("LEP FBA TEST")
# println()

# # EX_gal_e
# let 
#     _core_net0 = deepcopy(core_net0)
#     biom_id = extras(_core_net0, "BIOM")
#     linear_weights!(_core_net0, biom_id, 1.0)
#     exchs = ["EX_glc__D_e", "EX_lac__D_e", "EX_malt_e", "EX_gal_e", "EX_glyc_e", "EX_ac_e"]
#     lb!(_core_net0, exchs, 0.0)
    
#     for nut_id in exchs

#         println("-"^40)
#         println(nut_id)

#         lb!(_core_net0, nut_id, -10.0)
#         opm = fba(_core_net0, Clp.Optimizer)
#         println(biom_id, ": ", solution(opm, biom_id))
#         # println(nut_id, ": ", solution(opm, nut_id))
#         @assert solution(opm, biom_id) > 1e-2
#         lb!(_core_net0, nut_id, 0.0)
#     end    

#     # -------------------------------------------
#     # EXP

#     println()
#     println("="^40)
#     println("EXPERIMENTS")
#     println()
    
#     # -------------------------------------------
#     println("-"^40)
#     println("Glc phase")
#     lb!(_core_net0, exchs, 0.0)
#     exch = "EX_glc__D_e"
#     lb!(_core_net0, exch, lb(core_net0, exch))
#     opm = fba(_core_net0, Clp.Optimizer)
#     println(biom_id, ": ", solution(opm, biom_id))
#     lb!(_core_net0, exchs, 0.0)
    
#     # -------------------------------------------
#     println("-"^40)
#     println("Lac-Mal-Gal phase")
#     lb!(_core_net0, exchs, 0.0)
#     exch = "EX_lac__D_e"
#     lb!(_core_net0, exch, lb(core_net0, exch))
#     exch = "EX_malt_e"
#     lb!(_core_net0, exch, lb(core_net0, exch))
#     exch = "EX_gal_e"
#     lb!(_core_net0, exch, lb(core_net0, exch))
#     opm = fba(_core_net0, Clp.Optimizer)
#     println(biom_id, ": ", solution(opm, biom_id))
#     lb!(_core_net0, exchs, 0.0)

#     # -------------------------------------------
#     println("-"^40)
#     println("Glyc-Ac phase")
#     lb!(_core_net0, exchs, 0.0)
#     exch = "EX_glyc_e"
#     lb!(_core_net0, exch, lb(core_net0, exch))
#     exch = "EX_ac_e"
#     lb!(_core_net0, exch, lb(core_net0, exch))
#     opm = fba(_core_net0, Clp.Optimizer)
#     println(biom_id, ": ", solution(opm, biom_id))
#     lb!(_core_net0, exchs, 0.0)
# end
