## ------------------------------------------------------------
@time begin
    using MetXGEMs
    using MetXBase
    using MetXOptim
    using MetXNetHub
    using NutrientLimitedGEMs
end

# ------------------------------------------------------------
include("1_setup.jl")
include("2_utils.jl")

## ------------------------------------------------------------
# Prepare network
@tempcontext ["GEM_CORE_RXN_MAP" => v"0.1.0"] let
    @stage! rxn_map = Dict()
    
    # [1] ACALD::[Acetaldehyde dehydrogenase (acetylating)] Acetaldehyde dehydrogenase (acetylating)
    # (-1.0) nad_c::[Nicotinamide adenine dinucleotide] + (-1.0) acald_c::[Acetaldehyde] + (-1.0) coa_c::[Coenzyme A] <==> (1.0) h_c::[H+] + (1.0) nadh_c::[Nicotinamide adenine dinucleotide - reduced] + (1.0) accoa_c::[Acetyl-CoA]
    rxn_map["ACALD"] = "ACALD"

    # [2] ACALDt::[Acetaldehyde reversible transport] Acetaldehyde reversible transport
    # (-1.0) acald_e::[Acetaldehyde] <==> (1.0) acald_c::[Acetaldehyde]
    rxn_map["ACALDt"] = "ACALDtpp"

    # [3] ACKr::[Acetate kinase] Acetate kinase
    # (-1.0) ac_c::[Acetate] + (-1.0) atp_c::[ATP C10H12N5O13P3] <==> (1.0) actp_c::[Acetyl phosphate] + (1.0) adp_c::[ADP C10H12N5O10P2]
    rxn_map["ACKr"] = "ACKr"

    # [4] ACONTa::[Aconitase (half-reaction A, Citrate hydro-lyase)] Aconitase (half-reaction A, Citrate hydro-lyase)
    # (-1.0) cit_c::[Citrate] <==> (1.0) h2o_c::[H2O H2O] + (1.0) acon_C_c::[Cis-Aconitate]
    rxn_map["ACONTa"] = "ACONTa"

    # [5] ACONTb::[Aconitase (half-reaction B, Isocitrate hydro-lyase)] Aconitase (half-reaction B, Isocitrate hydro-lyase)
    # (-1.0) h2o_c::[H2O H2O] + (-1.0) acon_C_c::[Cis-Aconitate] <==> (1.0) icit_c::[Isocitrate]
    rxn_map["ACONTb"] = "ACONTb"

    # [6] ACt2r::[Acetate reversible transport via proton symport] Acetate reversible transport via proton symport
    # (-1.0) h_e::[H+] + (-1.0) ac_e::[Acetate] <==> (1.0) h_c::[H+] + (1.0) ac_c::[Acetate]
    rxn_map["ACt2r"] = "ACt2rpp"

    # [7] ADK1::[Adenylate kinase] Adenylate kinase
    # (-1.0) amp_c::[AMP C10H12N5O7P] + (-1.0) atp_c::[ATP C10H12N5O13P3] <==> (2.0) adp_c::[ADP C10H12N5O10P2]
    rxn_map["ADK1"] = "ADK1"

    # [8] AKGDH::[2-Oxogluterate dehydrogenase] 2-Oxogluterate dehydrogenase
    # (-1.0) nad_c::[Nicotinamide adenine dinucleotide] + (-1.0) akg_c::[2-Oxoglutarate] + (-1.0) coa_c::[Coenzyme A] ==> (1.0) nadh_c::[Nicotinamide adenine dinucleotide - reduced] + (1.0) succoa_c::[Succinyl-CoA] + (1.0) co2_c::[CO2 CO2]
    rxn_map["AKGDH"] = "AKGDH"

    # [9] AKGt2r::[2 oxoglutarate reversible transport via symport] 2 oxoglutarate reversible transport via symport
    # (-1.0) h_e::[H+] + (-1.0) akg_e::[2-Oxoglutarate] <==> (1.0) h_c::[H+] + (1.0) akg_c::[2-Oxoglutarate]
    rxn_map["AKGt2r"] = "AKGt2rpp"

    # [10] ALCD2x::[Alcohol dehydrogenase (ethanol)] Alcohol dehydrogenase (ethanol)
    # (-1.0) nad_c::[Nicotinamide adenine dinucleotide] + (-1.0) etoh_c::[Ethanol] <==> (1.0) h_c::[H+] + (1.0) nadh_c::[Nicotinamide adenine dinucleotide - reduced] + (1.0) acald_c::[Acetaldehyde]
    rxn_map["ALCD2x"] = "ALCD2x"

    # [11] AMALT1::[Amylomaltase (maltotriose)] Alternate Carbon Metabolism
    # (-1.0) malt_c::[Maltose C12H22O11] ==> (2.0) glc__D_e::[D-Glucose]
    rxn_map["AMALT1"] = "AMALT1"

    # [12] ATPM::[ATP maintenance requirement] ATP maintenance requirement
    # (-1.0) h2o_c::[H2O H2O] + (-1.0) atp_c::[ATP C10H12N5O13P3] ==> (1.0) h_c::[H+] + (1.0) pi_c::[Phosphate] + (1.0) adp_c::[ADP C10H12N5O10P2]
    rxn_map["ATPM"] = "ATPM"

    # [13] ATPS4r::[ATP synthase (four protons for one ATP)] ATP synthase (four protons for one ATP)
    # (-4.0) h_e::[H+] + (-1.0) pi_c::[Phosphate] + (-1.0) adp_c::[ADP C10H12N5O10P2] <==> (1.0) h2o_c::[H2O H2O] + (3.0) h_c::[H+] + (1.0) atp_c::[ATP C10H12N5O13P3]
    rxn_map["ATPS4r"] = "ATPS4rpp"

    # [14] BIOMASS_Ecoli_core_w_GAM::[Biomass Objective Function with GAM] Biomass Objective Function with GAM
    # (-0.2557) gln__L_c::[L-Glutamine] + (-4.9414) glu__L_c::[L-Glutamate] + (-59.81) h2o_c::[H2O H2O] + (-3.547) nad_c::[Nicotinamide adenine dinucleotide] + (-13.0279) nadph_c::[Nicotinamide adenine dinucleotide phosphate - reduced] + (-1.496) 3pg_c::[3-Phospho-D-glycerate] + (-1.7867) oaa_c::[Oxaloacetate] + (-0.5191) pep_c::[Phosphoenolpyruvate] + (-2.8328) pyr_c::[Pyruvate] + (-0.8977) r5p_c::[Alpha-D-Ribose 5-phosphate] + (-3.7478) accoa_c::[Acetyl-CoA] + (-59.81) atp_c::[ATP C10H12N5O13P3] + (-0.361) e4p_c::[D-Erythrose 4-phosphate] + (-0.0709) f6p_c::[D-Fructose 6-phosphate] + (-0.129) g3p_c::[Glyceraldehyde 3-phosphate] + (-0.205) g6p_c::[D-Glucose 6-phosphate] ==> (59.81) h_c::[H+] + (3.547) nadh_c::[Nicotinamide adenine dinucleotide - reduced] + (13.0279) nadp_c::[Nicotinamide adenine dinucleotide phosphate] + (59.81) pi_c::[Phosphate] + (59.81) adp_c::[ADP C10H12N5O10P2] + (4.1182) akg_c::[2-Oxoglutarate] + (3.7478) coa_c::[Coenzyme A]
    rxn_map["BIOMASS_Ecoli_core_w_GAM"] = "BIOMASS_Ec_iJO1366_core_53p95M"

    # [15] CO2t::[CO2 transporter via diffusion] CO2 transporter via diffusion
    # (-1.0) co2_e::[CO2 CO2] <==> (1.0) co2_c::[CO2 CO2]
    rxn_map["CO2t"] = "CO2tpp"

    # [16] CS::[Citrate synthase] Citrate synthase
    # (-1.0) h2o_c::[H2O H2O] + (-1.0) oaa_c::[Oxaloacetate] + (-1.0) accoa_c::[Acetyl-CoA] ==> (1.0) h_c::[H+] + (1.0) cit_c::[Citrate] + (1.0) coa_c::[Coenzyme A]
    rxn_map["CS"] = "CS"

    # [17] CYTBD::[Cytochrome oxidase bd (ubiquinol-8: 2 protons)] Cytochrome oxidase bd (ubiquinol-8: 2 protons)
    # (-2.0) h_c::[H+] + (-0.5) o2_c::[O2 O2] + (-1.0) q8h2_c::[Ubiquinol-8] ==> (1.0) h2o_c::[H2O H2O] + (2.0) h_e::[H+] + (1.0) q8_c::[Ubiquinone-8]
    rxn_map["CYTBD"] = "CYTBDpp"

    # [18] D_LACt2::[D lactate transport via proton symport] D lactate transport via proton symport
    # (-1.0) h_e::[H+] + (-1.0) lac__D_e::[D-Lactate] <==> (1.0) h_c::[H+] + (1.0) lac__D_c::[D-Lactate]
    rxn_map["D_LACt2"] = "D_LACt2pp"

    # [19] ENO::[Enolase] Enolase
    # (-1.0) 2pg_c::[D-Glycerate 2-phosphate] <==> (1.0) h2o_c::[H2O H2O] + (1.0) pep_c::[Phosphoenolpyruvate]
    rxn_map["ENO"] = "ENO"

    # [20] ETOHt2r::[Ethanol reversible transport via proton symport] Ethanol reversible transport via proton symport
    # (-1.0) h_e::[H+] + (-1.0) etoh_e::[Ethanol] <==> (1.0) h_c::[H+] + (1.0) etoh_c::[Ethanol]
    rxn_map["ETOHt2r"] = "ETOHtrpp"

    # [21] EX_ac_e::[Acetate exchange] Acetate exchange
    # (-1.0) ac_e::[Acetate] <==> 
    rxn_map["EX_ac_e"] = "EX_ac_e"

    # [22] EX_acald_e::[Acetaldehyde exchange] Acetaldehyde exchange
    # (-1.0) acald_e::[Acetaldehyde] ==> 
    rxn_map["EX_acald_e"] = "EX_acald_e"

    # [23] EX_akg_e::[2-Oxoglutarate exchange] 2-Oxoglutarate exchange
    # (-1.0) akg_e::[2-Oxoglutarate] ==> 
    rxn_map["EX_akg_e"] = "EX_akg_e"

    # [24] EX_co2_e::[CO2 exchange] CO2 exchange
    # (-1.0) co2_e::[CO2 CO2] <==> 
    rxn_map["EX_co2_e"] = "EX_co2_e"

    # [25] EX_etoh_e::[Ethanol exchange] Ethanol exchange
    # (-1.0) etoh_e::[Ethanol] ==> 
    rxn_map["EX_etoh_e"] = "EX_etoh_e"

    # [26] EX_for_e::[Formate exchange] Formate exchange
    # (-1.0) for_e::[Formate] ==> 
    rxn_map["EX_for_e"] = "EX_for_e"

    # [27] EX_fru_e::[D-Fructose exchange] D-Fructose exchange
    # (-1.0) fru_e::[D-Fructose] ==> 
    rxn_map["EX_fru_e"] = "EX_fru_e"

    # [28] EX_fum_e::[Fumarate exchange] Fumarate exchange
    # (-1.0) fum_e::[Fumarate] ==> 
    rxn_map["EX_fum_e"] = "EX_fum_e"

    # [29] EX_gal_e::[D-Galactose exchange] Extracellular exchange
    # (-1.0) gal_e::[D-Galactose] <==> 
    rxn_map["EX_gal_e"] = "EX_gal_e"

    # [30] EX_glc__D_e::[D-Glucose exchange] D-Glucose exchange
    # (-1.0) glc__D_e::[D-Glucose] <==> 
    rxn_map["EX_glc__D_e"] = "EX_glc__D_e"

    # [31] EX_gln__L_e::[L-Glutamine exchange] L-Glutamine exchange
    # (-1.0) gln__L_e::[L-Glutamine] ==> 
    rxn_map["EX_gln__L_e"] = "EX_gln__L_e"

    # [32] EX_glu__L_e::[L-Glutamate exchange] L-Glutamate exchange
    # (-1.0) glu__L_e::[L-Glutamate] ==> 
    rxn_map["EX_glu__L_e"] = "EX_glu__L_e"

    # [33] EX_glyc_e::[Glycerol exchange] Extracellular exchange
    # (-1.0) glyc_e::[Glycerol] <==> 
    rxn_map["EX_glyc_e"] = "EX_glyc_e"

    # [34] EX_h2o_e::[H2O exchange] H2O exchange
    # (-1.0) h2o_e::[H2O H2O] <==> 
    rxn_map["EX_h2o_e"] = "EX_h2o_e"

    # [35] EX_h_e::[H+ exchange] H+ exchange
    # (-1.0) h_e::[H+] <==> 
    rxn_map["EX_h_e"] = "EX_h_e"

    # [36] EX_lac__D_e::[D-lactate exchange] D-lactate exchange
    # (-1.0) lac__D_e::[D-Lactate] <==> 
    rxn_map["EX_lac__D_e"] = "EX_lac__D_e"

    # [37] EX_mal__L_e::[L-Malate exchange] L-Malate exchange
    # (-1.0) mal__L_e::[L-Malate] ==> 
    rxn_map["EX_mal__L_e"] = "EX_mal__L_e"

    # [38] EX_malt_e::[Maltose exchange] Extracellular exchange
    # (-1.0) malt_e::[Maltose C12H22O11] <==> 
    rxn_map["EX_malt_e"] = "EX_malt_e"

    # [39] EX_nh4_e::[Ammonia exchange] Ammonia exchange
    # (-1.0) nh4_e::[Ammonium] <==> 
    rxn_map["EX_nh4_e"] = "EX_nh4_e"

    # [40] EX_o2_e::[O2 exchange] O2 exchange
    # (-1.0) o2_e::[O2 O2] <==> 
    rxn_map["EX_o2_e"] = "EX_o2_e"

    # [41] EX_pi_e::[Phosphate exchange] Phosphate exchange
    # (-1.0) pi_e::[Phosphate] <==> 
    rxn_map["EX_pi_e"] = "EX_pi_e"

    # [42] EX_pyr_e::[Pyruvate exchange] Pyruvate exchange
    # (-1.0) pyr_e::[Pyruvate] ==> 
    rxn_map["EX_pyr_e"] = "EX_pyr_e"

    # [43] EX_succ_e::[Succinate exchange] Succinate exchange
    # (-1.0) succ_e::[Succinate] ==> 
    rxn_map["EX_succ_e"] = "EX_succ_e"

    # [44] FBA::[Fructose-bisphosphate aldolase] Fructose-bisphosphate aldolase
    # (-1.0) fdp_c::[D-Fructose 1,6-bisphosphate] <==> (1.0) dhap_c::[Dihydroxyacetone phosphate] + (1.0) g3p_c::[Glyceraldehyde 3-phosphate]
    rxn_map["FBA"] = "FBA"

    # [45] FBP::[Fructose-bisphosphatase] Fructose-bisphosphatase
    # (-1.0) h2o_c::[H2O H2O] + (-1.0) fdp_c::[D-Fructose 1,6-bisphosphate] ==> (1.0) pi_c::[Phosphate] + (1.0) f6p_c::[D-Fructose 6-phosphate]
    rxn_map["FBP"] = "FBP"

    # [46] FORt::[Formate transport via diffusion] Formate transport via diffusion
    # (-1.0) for_e::[Formate] <== (1.0) for_c::[Formate]
    rxn_map["FORt"] = "FORtppi"

    # [47] FORt2::[Formate transport in via proton symport] Formate transport in via proton symport
    # (-1.0) h_e::[H+] + (-1.0) for_e::[Formate] ==> (1.0) h_c::[H+] + (1.0) for_c::[Formate]
    rxn_map["FORt2"] = "FORt2pp"

    # [48] FRD7::[Fumarate reductase] Fumarate reductase
    # (-1.0) q8h2_c::[Ubiquinol-8] + (-1.0) fum_c::[Fumarate] ==> (1.0) q8_c::[Ubiquinone-8] + (1.0) succ_c::[Succinate]
    rxn_map["FRD7"] = "FRD2"

    # [49] FRUpts2::[Fructose transport via PEP:Pyr PTS (f6p generating)] Fructose transport via PEP:Pyr PTS (f6p generating)
    # (-1.0) pep_c::[Phosphoenolpyruvate] + (-1.0) fru_e::[D-Fructose] ==> (1.0) pyr_c::[Pyruvate] + (1.0) f6p_c::[D-Fructose 6-phosphate]
    rxn_map["FRUpts2"] = "FRUpts2pp"

    # [50] FUM::[Fumarase] Fumarase
    # (-1.0) h2o_c::[H2O H2O] + (-1.0) fum_c::[Fumarate] <==> (1.0) mal__L_c::[L-Malate]
    rxn_map["FUM"] = "FUM"

    # [51] FUMt2_2::[Fumarate transport via proton symport (2 H)] Fumarate transport via proton symport (2 H)
    # (-2.0) h_e::[H+] + (-1.0) fum_e::[Fumarate] ==> (2.0) h_c::[H+] + (1.0) fum_c::[Fumarate]
    rxn_map["FUMt2_2"] = "FUMt2_2pp"

    # [52] G3PD2::[Glycerol-3-phosphate dehydrogenase (NADP)] Alternate Carbon Metabolism
    # (-1.0) nadp_c::[Nicotinamide adenine dinucleotide phosphate] + (-1.0) glyc3p_c::[Glycerol 3-phosphate] <==> (1.0) h_c::[H+] + (1.0) nadph_c::[Nicotinamide adenine dinucleotide phosphate - reduced] + (1.0) dhap_c::[Dihydroxyacetone phosphate]
    rxn_map["G3PD2"] = "G3PD2"

    # [53] G6PDH2r::[Glucose 6-phosphate dehydrogenase] Glucose 6-phosphate dehydrogenase
    # (-1.0) nadp_c::[Nicotinamide adenine dinucleotide phosphate] + (-1.0) g6p_c::[D-Glucose 6-phosphate] <==> (1.0) h_c::[H+] + (1.0) nadph_c::[Nicotinamide adenine dinucleotide phosphate - reduced] + (1.0) 6pgl_c::[6-phospho-D-glucono-1,5-lactone]
    rxn_map["G6PDH2r"] = "G6PDH2r"

    # [54] GALKr::[Galactokinase] Alternate Carbon Metabolism
    # (-1.0) atp_c::[ATP C10H12N5O13P3] + (-1.0) gal_c::[D-Galactose] <==> (1.0) h_c::[H+] + (1.0) adp_c::[ADP C10H12N5O10P2] + (1.0) gal1p_c::[Alpha-D-Galactose 1-phosphate]
    rxn_map["GALKr"] = "GALKr"

    # [55] GALabcpp::[D-galactose transport via ABC system (periplasm)] Transport, Inner Membrane
    # (-1.0) h2o_e::[H2O H2O] + (-1.0) atp_c::[ATP C10H12N5O13P3] + (-1.0) gal_e::[D-Galactose] ==> (1.0) h_c::[H+] + (1.0) pi_c::[Phosphate] + (1.0) adp_c::[ADP C10H12N5O10P2] + (1.0) gal_c::[D-Galactose]
    rxn_map["GALabcpp"] = "GALabcpp"

    # [56] GAPD::[Glyceraldehyde-3-phosphate dehydrogenase] Glyceraldehyde-3-phosphate dehydrogenase
    # (-1.0) nad_c::[Nicotinamide adenine dinucleotide] + (-1.0) pi_c::[Phosphate] + (-1.0) g3p_c::[Glyceraldehyde 3-phosphate] <==> (1.0) h_c::[H+] + (1.0) nadh_c::[Nicotinamide adenine dinucleotide - reduced] + (1.0) 13dpg_c::[3-Phospho-D-glyceroyl phosphate]
    rxn_map["GAPD"] = "GAPD"

    # [57] GLCpts::[D-glucose transport via PEP:Pyr PTS] D-glucose transport via PEP:Pyr PTS
    # (-1.0) glc__D_e::[D-Glucose] + (-1.0) pep_c::[Phosphoenolpyruvate] ==> (1.0) pyr_c::[Pyruvate] + (1.0) g6p_c::[D-Glucose 6-phosphate]
    rxn_map["GLCpts"] = "GLCptspp"

    # [58] GLNS::[Glutamine synthetase] Glutamine synthetase
    # (-1.0) glu__L_c::[L-Glutamate] + (-1.0) nh4_c::[Ammonium] + (-1.0) atp_c::[ATP C10H12N5O13P3] ==> (1.0) gln__L_c::[L-Glutamine] + (1.0) h_c::[H+] + (1.0) pi_c::[Phosphate] + (1.0) adp_c::[ADP C10H12N5O10P2]
    rxn_map["GLNS"] = "GLNS"

    # [59] GLNabc::[L-glutamine transport via ABC system] L-glutamine transport via ABC system
    # (-1.0) gln__L_e::[L-Glutamine] + (-1.0) h2o_c::[H2O H2O] + (-1.0) atp_c::[ATP C10H12N5O13P3] ==> (1.0) gln__L_c::[L-Glutamine] + (1.0) h_c::[H+] + (1.0) pi_c::[Phosphate] + (1.0) adp_c::[ADP C10H12N5O10P2]
    rxn_map["GLNabc"] = "GLNabcpp"

    # [60] GLUDy::[Glutamate dehydrogenase (NADP)] Glutamate dehydrogenase (NADP)
    # (-1.0) glu__L_c::[L-Glutamate] + (-1.0) h2o_c::[H2O H2O] + (-1.0) nadp_c::[Nicotinamide adenine dinucleotide phosphate] <==> (1.0) h_c::[H+] + (1.0) nadph_c::[Nicotinamide adenine dinucleotide phosphate - reduced] + (1.0) nh4_c::[Ammonium] + (1.0) akg_c::[2-Oxoglutarate]
    rxn_map["GLUDy"] = "GLUDy"

    # [61] GLUN::[Glutaminase] Glutaminase
    # (-1.0) gln__L_c::[L-Glutamine] + (-1.0) h2o_c::[H2O H2O] ==> (1.0) glu__L_c::[L-Glutamate] + (1.0) nh4_c::[Ammonium]
    rxn_map["GLUN"] = "GLUN"

    # [62] GLUSy::[Glutamate synthase (NADPH)] Glutamate synthase (NADPH)
    # (-1.0) gln__L_c::[L-Glutamine] + (-1.0) h_c::[H+] + (-1.0) nadph_c::[Nicotinamide adenine dinucleotide phosphate - reduced] + (-1.0) akg_c::[2-Oxoglutarate] ==> (2.0) glu__L_c::[L-Glutamate] + (1.0) nadp_c::[Nicotinamide adenine dinucleotide phosphate]
    rxn_map["GLUSy"] = "GLUSy"

    # [63] GLUt2r::[L glutamate transport via proton symport  reversible] L glutamate transport via proton symport  reversible
    # (-1.0) glu__L_e::[L-Glutamate] + (-1.0) h_e::[H+] <==> (1.0) glu__L_c::[L-Glutamate] + (1.0) h_c::[H+]
    rxn_map["GLUt2r"] = "GLUt2rpp"

    # [64] GLYCtpp::[Glycerol transport via channel (periplasm)] Transport, Inner Membrane
    # (-1.0) glyc_e::[Glycerol] <==> (1.0) glyc_c::[Glycerol]
    rxn_map["GLYCtpp"] = "GLYCtpp"

    # [65] GLYK::[Glycerol kinase] Alternate Carbon Metabolism
    # (-1.0) atp_c::[ATP C10H12N5O13P3] + (-1.0) glyc_c::[Glycerol] ==> (1.0) h_c::[H+] + (1.0) adp_c::[ADP C10H12N5O10P2] + (1.0) glyc3p_c::[Glycerol 3-phosphate]
    rxn_map["GLYK"] = "GLYK"

    # [66] GND::[Phosphogluconate dehydrogenase] Phosphogluconate dehydrogenase
    # (-1.0) nadp_c::[Nicotinamide adenine dinucleotide phosphate] + (-1.0) 6pgc_c::[6-Phospho-D-gluconate] ==> (1.0) nadph_c::[Nicotinamide adenine dinucleotide phosphate - reduced] + (1.0) ru5p__D_c::[D-Ribulose 5-phosphate] + (1.0) co2_c::[CO2 CO2]
    rxn_map["GND"] = "GND"

    # [67] H2Ot::[H2O transport via diffusion] H2O transport via diffusion
    # (-1.0) h2o_e::[H2O H2O] <==> (1.0) h2o_c::[H2O H2O]
    rxn_map["H2Ot"] = "H2Otex"

    # [68] ICDHyr::[Isocitrate dehydrogenase (NADP)] Isocitrate dehydrogenase (NADP)
    # (-1.0) icit_c::[Isocitrate] + (-1.0) nadp_c::[Nicotinamide adenine dinucleotide phosphate] <==> (1.0) nadph_c::[Nicotinamide adenine dinucleotide phosphate - reduced] + (1.0) akg_c::[2-Oxoglutarate] + (1.0) co2_c::[CO2 CO2]
    rxn_map["ICDHyr"] = "ICDHyr"

    # [69] ICL::[Isocitrate lyase] Isocitrate lyase
    # (-1.0) icit_c::[Isocitrate] ==> (1.0) glx_c::[Glyoxylate] + (1.0) succ_c::[Succinate]
    rxn_map["ICL"] = "ICL"

    # [70] LDH_D::[D-lactate dehydrogenase] D-lactate dehydrogenase
    # (-1.0) lac__D_c::[D-Lactate] + (-1.0) nad_c::[Nicotinamide adenine dinucleotide] <==> (1.0) h_c::[H+] + (1.0) nadh_c::[Nicotinamide adenine dinucleotide - reduced] + (1.0) pyr_c::[Pyruvate]
    rxn_map["LDH_D"] = "LDH_D"

    # [71] MALS::[Malate synthase] Malate synthase
    # (-1.0) glx_c::[Glyoxylate] + (-1.0) h2o_c::[H2O H2O] + (-1.0) accoa_c::[Acetyl-CoA] ==> (1.0) h_c::[H+] + (1.0) mal__L_c::[L-Malate] + (1.0) coa_c::[Coenzyme A]
    rxn_map["MALS"] = "MALS"

    # [72] MALTabcpp::[Maltose transport via ABC system (periplasm)] Transport, Inner Membrane
    # (-1.0) h2o_c::[H2O H2O] + (-1.0) atp_c::[ATP C10H12N5O13P3] + (-1.0) malt_e::[Maltose C12H22O11] ==> (1.0) h_c::[H+] + (1.0) pi_c::[Phosphate] + (1.0) adp_c::[ADP C10H12N5O10P2] + (1.0) malt_c::[Maltose C12H22O11]
    rxn_map["MALTabcpp"] = "MALTabcpp"

    # [73] MALt2_2::[Malate transport via proton symport (2 H)] Malate transport via proton symport (2 H)
    # (-2.0) h_e::[H+] + (-1.0) mal__L_e::[L-Malate] ==> (2.0) h_c::[H+] + (1.0) mal__L_c::[L-Malate]
    rxn_map["MALt2_2"] = "MALt2_2pp"

    # [74] MDH::[Malate dehydrogenase] Malate dehydrogenase
    # (-1.0) mal__L_c::[L-Malate] + (-1.0) nad_c::[Nicotinamide adenine dinucleotide] <==> (1.0) h_c::[H+] + (1.0) nadh_c::[Nicotinamide adenine dinucleotide - reduced] + (1.0) oaa_c::[Oxaloacetate]
    rxn_map["MDH"] = "MDH"

    # [75] ME1::[Malic enzyme (NAD)] Malic enzyme (NAD)
    # (-1.0) mal__L_c::[L-Malate] + (-1.0) nad_c::[Nicotinamide adenine dinucleotide] ==> (1.0) nadh_c::[Nicotinamide adenine dinucleotide - reduced] + (1.0) pyr_c::[Pyruvate] + (1.0) co2_c::[CO2 CO2]
    rxn_map["ME1"] = "ME1"

    # [76] ME2::[Malic enzyme (NADP)] Malic enzyme (NADP)
    # (-1.0) mal__L_c::[L-Malate] + (-1.0) nadp_c::[Nicotinamide adenine dinucleotide phosphate] ==> (1.0) nadph_c::[Nicotinamide adenine dinucleotide phosphate - reduced] + (1.0) pyr_c::[Pyruvate] + (1.0) co2_c::[CO2 CO2]
    rxn_map["ME2"] = "ME2"

    # [77] NADH16::[NADH dehydrogenase (ubiquinone-8 & 3 protons)] NADH dehydrogenase (ubiquinone-8 & 3 protons)
    # (-4.0) h_c::[H+] + (-1.0) nadh_c::[Nicotinamide adenine dinucleotide - reduced] + (-1.0) q8_c::[Ubiquinone-8] ==> (3.0) h_e::[H+] + (1.0) nad_c::[Nicotinamide adenine dinucleotide] + (1.0) q8h2_c::[Ubiquinol-8]
    rxn_map["NADH16"] = "NADH16pp"

    # [78] NADTRHD::[NAD transhydrogenase] NAD transhydrogenase
    # (-1.0) nad_c::[Nicotinamide adenine dinucleotide] + (-1.0) nadph_c::[Nicotinamide adenine dinucleotide phosphate - reduced] ==> (1.0) nadh_c::[Nicotinamide adenine dinucleotide - reduced] + (1.0) nadp_c::[Nicotinamide adenine dinucleotide phosphate]
    rxn_map["NADTRHD"] = "NADTRHD"

    # [79] NH4t::[Ammonia reversible transport] Ammonia reversible transport
    # (-1.0) nh4_e::[Ammonium] <==> (1.0) nh4_c::[Ammonium]
    rxn_map["NH4t"] = "NH4tpp"

    # [80] O2t::[O2 transport diffusion ] O2 transport  diffusion 
    # (-1.0) o2_e::[O2 O2] <==> (1.0) o2_c::[O2 O2]
    rxn_map["O2t"] = "O2tpp"

    # [81] PDH::[Pyruvate dehydrogenase] Pyruvate dehydrogenase
    # (-1.0) nad_c::[Nicotinamide adenine dinucleotide] + (-1.0) pyr_c::[Pyruvate] + (-1.0) coa_c::[Coenzyme A] ==> (1.0) nadh_c::[Nicotinamide adenine dinucleotide - reduced] + (1.0) accoa_c::[Acetyl-CoA] + (1.0) co2_c::[CO2 CO2]
    rxn_map["PDH"] = "PDH"

    # [82] PFK::[Phosphofructokinase] Phosphofructokinase
    # (-1.0) atp_c::[ATP C10H12N5O13P3] + (-1.0) f6p_c::[D-Fructose 6-phosphate] ==> (1.0) h_c::[H+] + (1.0) adp_c::[ADP C10H12N5O10P2] + (1.0) fdp_c::[D-Fructose 1,6-bisphosphate]
    rxn_map["PFK"] = "PFK"

    # [83] PFL::[Pyruvate formate lyase] Pyruvate formate lyase
    # (-1.0) pyr_c::[Pyruvate] + (-1.0) coa_c::[Coenzyme A] ==> (1.0) accoa_c::[Acetyl-CoA] + (1.0) for_c::[Formate]
    rxn_map["PFL"] = 

    # [84] PGI::[Glucose-6-phosphate isomerase] Glucose-6-phosphate isomerase
    # (-1.0) g6p_c::[D-Glucose 6-phosphate] <==> (1.0) f6p_c::[D-Fructose 6-phosphate]
    rxn_map["PGI"] = "PGI"

    # [85] PGK::[Phosphoglycerate kinase] Phosphoglycerate kinase
    # (-1.0) 3pg_c::[3-Phospho-D-glycerate] + (-1.0) atp_c::[ATP C10H12N5O13P3] <==> (1.0) 13dpg_c::[3-Phospho-D-glyceroyl phosphate] + (1.0) adp_c::[ADP C10H12N5O10P2]
    rxn_map["PGK"] = "PGK"

    # [86] PGL::[6-phosphogluconolactonase] 6-phosphogluconolactonase
    # (-1.0) h2o_c::[H2O H2O] + (-1.0) 6pgl_c::[6-phospho-D-glucono-1,5-lactone] ==> (1.0) h_c::[H+] + (1.0) 6pgc_c::[6-Phospho-D-gluconate]
    rxn_map["PGL"] = "PGL"

    # [87] PGM::[Phosphoglycerate mutase] Phosphoglycerate mutase
    # (-1.0) 2pg_c::[D-Glycerate 2-phosphate] <==> (1.0) 3pg_c::[3-Phospho-D-glycerate]
    rxn_map["PGM"] = "PGM"

    # [88] PGMT::[Phosphoglucomutase] Alternate Carbon Metabolism
    # (-1.0) g1p_c::[D-Glucose 1-phosphate] <==> (1.0) g6p_c::[D-Glucose 6-phosphate]
    rxn_map["PGMT"] = "PGMT"

    # [89] PIt2r::[Phosphate reversible transport via symport] Phosphate reversible transport via symport
    # (-1.0) h_e::[H+] + (-1.0) pi_e::[Phosphate] <==> (1.0) h_c::[H+] + (1.0) pi_c::[Phosphate]
    rxn_map["PIt2r"] = "PIt2rpp"

    # [90] PPC::[Phosphoenolpyruvate carboxylase] Phosphoenolpyruvate carboxylase
    # (-1.0) h2o_c::[H2O H2O] + (-1.0) pep_c::[Phosphoenolpyruvate] + (-1.0) co2_c::[CO2 CO2] ==> (1.0) h_c::[H+] + (1.0) oaa_c::[Oxaloacetate] + (1.0) pi_c::[Phosphate]
    rxn_map["PPC"] = "PPC"

    # [91] PPCK::[Phosphoenolpyruvate carboxykinase] Phosphoenolpyruvate carboxykinase
    # (-1.0) oaa_c::[Oxaloacetate] + (-1.0) atp_c::[ATP C10H12N5O13P3] ==> (1.0) pep_c::[Phosphoenolpyruvate] + (1.0) adp_c::[ADP C10H12N5O10P2] + (1.0) co2_c::[CO2 CO2]
    rxn_map["PPCK"] = "PPCK"

    # [92] PPS::[Phosphoenolpyruvate synthase] Phosphoenolpyruvate synthase
    # (-1.0) h2o_c::[H2O H2O] + (-1.0) pyr_c::[Pyruvate] + (-1.0) atp_c::[ATP C10H12N5O13P3] ==> (2.0) h_c::[H+] + (1.0) pep_c::[Phosphoenolpyruvate] + (1.0) pi_c::[Phosphate] + (1.0) amp_c::[AMP C10H12N5O7P]
    rxn_map["PPS"] = "PPS"

    # [93] PTAr::[Phosphotransacetylase] Phosphotransacetylase
    # (-1.0) pi_c::[Phosphate] + (-1.0) accoa_c::[Acetyl-CoA] <==> (1.0) actp_c::[Acetyl phosphate] + (1.0) coa_c::[Coenzyme A]
    rxn_map["PTAr"] = "PTAr"

    # [94] PYK::[Pyruvate kinase] Pyruvate kinase
    # (-1.0) h_c::[H+] + (-1.0) pep_c::[Phosphoenolpyruvate] + (-1.0) adp_c::[ADP C10H12N5O10P2] ==> (1.0) pyr_c::[Pyruvate] + (1.0) atp_c::[ATP C10H12N5O13P3]
    rxn_map["PYK"] = "PYK"

    # [95] PYRt2::[Pyruvate transport in via proton symport] Pyruvate transport in via proton symport
    # (-1.0) h_e::[H+] + (-1.0) pyr_e::[Pyruvate] <==> (1.0) h_c::[H+] + (1.0) pyr_c::[Pyruvate]
    rxn_map["PYRt2"] = "PYRtex"

    # [96] RPE::[Ribulose 5-phosphate 3-epimerase] Ribulose 5-phosphate 3-epimerase
    # (-1.0) ru5p__D_c::[D-Ribulose 5-phosphate] <==> (1.0) xu5p__D_c::[D-Xylulose 5-phosphate]
    rxn_map["RPE"] = "RPE"

    # [97] RPI::[Ribose-5-phosphate isomerase] Ribose-5-phosphate isomerase
    # (-1.0) r5p_c::[Alpha-D-Ribose 5-phosphate] <==> (1.0) ru5p__D_c::[D-Ribulose 5-phosphate]
    rxn_map["RPI"] = "RPI"

    # [98] SUCCt2_2::[Succinate transport via proton symport (2 H)] Succinate transport via proton symport (2 H)
    # (-2.0) h_e::[H+] + (-1.0) succ_e::[Succinate] ==> (2.0) h_c::[H+] + (1.0) succ_c::[Succinate]
    rxn_map["SUCCt2_2"] = "SUCCt2_2pp"

    # [99] SUCCt3::[Succinate transport out via proton antiport] Succinate transport out via proton antiport
    # (-1.0) h_e::[H+] + (-1.0) succ_c::[Succinate] ==> (1.0) h_c::[H+] + (1.0) succ_e::[Succinate]
    rxn_map["SUCCt3"] = "SUCCt3pp"

    # [100] SUCDi::[Succinate dehydrogenase (irreversible)] Succinate dehydrogenase (irreversible)
    # (-1.0) q8_c::[Ubiquinone-8] + (-1.0) succ_c::[Succinate] ==> (1.0) q8h2_c::[Ubiquinol-8] + (1.0) fum_c::[Fumarate]
    rxn_map["SUCDi"] = "SUCDi"

    # [101] SUCOAS::[Succinyl-CoA synthetase (ADP-forming)] Succinyl-CoA synthetase (ADP-forming)
    # (-1.0) succ_c::[Succinate] + (-1.0) atp_c::[ATP C10H12N5O13P3] + (-1.0) coa_c::[Coenzyme A] <==> (1.0) pi_c::[Phosphate] + (1.0) succoa_c::[Succinyl-CoA] + (1.0) adp_c::[ADP C10H12N5O10P2]
    rxn_map["SUCOAS"] = "SUCOAS"

    # [102] TALA::[Transaldolase] Transaldolase
    # (-1.0) s7p_c::[Sedoheptulose 7-phosphate] + (-1.0) g3p_c::[Glyceraldehyde 3-phosphate] <==> (1.0) e4p_c::[D-Erythrose 4-phosphate] + (1.0) f6p_c::[D-Fructose 6-phosphate]
    rxn_map["TALA"] = "TALA"

    # [103] THD2::[NAD(P) transhydrogenase] NAD(P) transhydrogenase
    # (-2.0) h_e::[H+] + (-1.0) nadh_c::[Nicotinamide adenine dinucleotide - reduced] + (-1.0) nadp_c::[Nicotinamide adenine dinucleotide phosphate] ==> (2.0) h_c::[H+] + (1.0) nad_c::[Nicotinamide adenine dinucleotide] + (1.0) nadph_c::[Nicotinamide adenine dinucleotide phosphate - reduced]
    rxn_map["THD2"] = "THD2pp"

    # [104] TKT1::[Transketolase] Transketolase
    # (-1.0) r5p_c::[Alpha-D-Ribose 5-phosphate] + (-1.0) xu5p__D_c::[D-Xylulose 5-phosphate] <==> (1.0) s7p_c::[Sedoheptulose 7-phosphate] + (1.0) g3p_c::[Glyceraldehyde 3-phosphate]
    rxn_map["TKT1"] = "TKT1"

    # [105] TKT2::[Transketolase] Transketolase
    # (-1.0) xu5p__D_c::[D-Xylulose 5-phosphate] + (-1.0) e4p_c::[D-Erythrose 4-phosphate] <==> (1.0) f6p_c::[D-Fructose 6-phosphate] + (1.0) g3p_c::[Glyceraldehyde 3-phosphate]
    rxn_map["TKT2"] = "TKT2"

    # [106] TPI::[Triose-phosphate isomerase] Triose-phosphate isomerase
    # (-1.0) dhap_c::[Dihydroxyacetone phosphate] <==> (1.0) g3p_c::[Glyceraldehyde 3-phosphate]
    rxn_map["TPI"] = "TPI"

    # [107] UDPG4E::[UDPglucose 4-epimerase] Alternate Carbon Metabolism
    # (-1.0) udpg_c::[UDPglucose] <==> (1.0) udpgal_c::[UDPgalactose]
    rxn_map["UDPG4E"] = "UDPG4E"

    # [108] UGLT::[UDPglucose--hexose-1-phosphate uridylyltransferase] Alternate Carbon Metabolism
    # (-1.0) gal1p_c::[Alpha-D-Galactose 1-phosphate] + (-1.0) udpg_c::[UDPglucose] <==> (1.0) udpgal_c::[UDPgalactose] + (1.0) g1p_c::[D-Glucose 1-phosphate]
    rxn_map["UGLT"] = "UGLT"

end

## ------------------------------------------------------------
# save
_save_contextdb(SIMVER)