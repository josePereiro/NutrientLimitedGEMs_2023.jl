## ------------------------------------------------------------
@time begin
    using MetXGEMs
    using MetXBase
    using MetXOptim
    using MetXNetHub
    using Statistics
    using NutrientLimitedGEMs
end

# ------------------------------------------------------------
include("1_setup_sim.jl")
include("1.1_utils.jl")

## ------------------------------------------------------------
# Prepare network
@tempcontext ["MICROARRAY_COMPLEX_MEDIUM" => v"0.1.0"] let

    @stage! raw_map = Dict()
    raw_map["Supp1"] = Dict()
    raw_map["Supp3"] = Dict()
    
    # [1] ACALD::[Acetaldehyde dehydrogenase (acetylating)] Acetaldehyde dehydrogenase (acetylating)
    # (-1.0) nad_c::[Nicotinamide adenine dinucleotide] + (-1.0) acald_c::[Acetaldehyde] + (-1.0) coa_c::[Coenzyme A] <==> (1.0) h_c::[H+] + (1.0) nadh_c::[Nicotinamide adenine dinucleotide - reduced] + (1.0) accoa_c::[Acetyl-CoA]
    raw_map["Supp1"]["ACALD"] = ["adhE"]
    raw_map["Supp3"]["ACALD"] = [
        Dict{String, Any}(
            "raw" => [5349.4, 5261.7, 6189.1, 6177.1, 5184.7, 8185.0, 8342.5, 8122.1, 6242.0, 4317.8, 10578.1, 12094.9], 
            "norm" => [-0.356688017, -0.380536108, -0.14633547, -0.149135416, -0.401804589, 0.256917315, 0.284414679, 0.245787688, -0.134056745, -0.665768683, 0.626943512, 0.820261834], 
            "desc" => "MG1655_adhE_b1241  DB_XREF=PID:g1787493  SEG=NC_000913:-1294669,1297344  LEN=2675  DEF=CoA-linked acetaldehyde dehydrogenase and iron-dependent alcohol dehydrogenase; pyruvate-formate-lyase deactivase",
        )
    ]

    # [2] ACALDt::[Acetaldehyde reversible transport] Acetaldehyde reversible transport
    # (-1.0) acald_e::[Acetaldehyde] <==> (1.0) acald_c::[Acetaldehyde]
    raw_map["Supp1"]["ACALDt"] = [] # Not found
    raw_map["Supp3"]["ACALDt"] = []


    # [3] ACKr::[Acetate kinase] Acetate kinase
    # (-1.0) ac_c::[Acetate] + (-1.0) atp_c::[ATP C10H12N5O13P3] <==> (1.0) actp_c::[Acetyl phosphate] + (1.0) adp_c::[ADP C10H12N5O10P2]
    raw_map["Supp1"]["ACKr"] = ["ackA"]
    raw_map["Supp3"]["ACKr"] = [
        Dict{String, Any}(
            "raw" => [3695.4, 3993.3, 4051.5, 4005.3, 7765.5, 2883.3, 2351.3, 4252.4, 3627.8, 3539.0, 5439.9, 4253.8], 
            "norm" => [-0.108018063, 0.003832863, 0.024707544, 0.008161709, 0.963330222, -0.466027644, -0.760289974, 0.094528713, -0.134653674, -0.170406835, 0.449831534, 0.095003607], 
            "desc" => "MG1655_ackA_b2296  DB_XREF=PID:g1788633  SEG=NC_000913:+2411490,2412692  LEN=1202  DEF=acetate kinase"
        ), 
    ]

    # [4] ACONTa::[Aconitase (half-reaction A, Citrate hydro-lyase)] Aconitase (half-reaction A, Citrate hydro-lyase)
    # (-1.0) cit_c::[Citrate] <==> (1.0) h2o_c::[H2O H2O] + (1.0) acon_C_c::[Cis-Aconitate]
    raw_map["Supp1"]["ACONTa"] = ["acnA"]
    raw_map["Supp3"]["ACONTa"] = [
        Dict{String, Any}(
            "raw" => [1769.1, 2219.6, 2195.1, 2468.4, 2222.9, 4210.3, 2720.0, 5770.8, 1529.8, 2412.1, 5687.3, 8090.4], 
            "norm" => [-0.761776896, -0.434492788, -0.450505831, -0.281216296, -0.432349448, 0.489130539, -0.141185844, 0.943978836, -0.971449443, -0.314502777, 0.922951411, 1.431418537], 
            "desc" => "MG1655_acnA_b1276  DB_XREF=PID:g1787531  SEG=NC_000913:+1333855,1336530  LEN=2675  DEF=aconitate hydrase 1"
        )
    ]

    # [5] ACONTb::[Aconitase (half-reaction B, Isocitrate hydro-lyase)] Aconitase (half-reaction B, Isocitrate hydro-lyase)
    # (-1.0) h2o_c::[H2O H2O] + (-1.0) acon_C_c::[Cis-Aconitate] <==> (1.0) icit_c::[Isocitrate]
    raw_map["Supp1"]["ACONTb"] = ["acnB"]
    raw_map["Supp3"]["ACONTb"] = [
        Dict{String, Any}(
            "raw" => [49.8, 52.6, 61.7, 81.1, 52.0, 38.0, 78.2, 90.1, 54.5, 89.3, 31.0, 67.3], 
            "norm" => [-0.249452493, -0.170535436, 0.059672254, 0.454103679, -0.187086612, -0.639598817, 0.401570372, 0.60592887, -0.119342006, 0.59306194, -0.93333002, 0.185008269], 
            "desc" => "CFT073_acnB_c0147  DB_XREF=26246058  SEG=NC_004431:+141792,144542  LEN=2750  DEF=Aconitate hydratase 2"
        ),
        Dict{String, Any}(
            "raw" => [6977.9, 5482.1, 4965.8, 3904.1, 4696.6, 5439.3, 3749.1, 1692.2, 1501.3, 1315.3, 1086.0, 629.0], 
            "norm" => [1.36056125, 1.012496972, 0.869794485, 0.522758334, 0.789385055, 1.001189326, 0.464312635, -0.683331583, -0.856019377, -1.046839778, -1.323207569, -2.11109975], 
            "desc" => "MG1655_acnB_b0118  DB_XREF=PID:g2367097  SEG=NC_000913:+131615,134212  LEN=2597  DEF=aconitate hydrase B"
        ),
        Dict{String, Any}(
            "raw" => [71.1, 21.5, 34.0, 33.4, 35.8, 56.3, 37.6, 22.5, 20.0, 13.8, 5.0, 1.3], 
            "norm" => [1.79863514, 0.07312224, 0.734320327, 0.708633683, 0.808745168, 1.461920503, 0.879518242, 0.138710582, -0.03121442, -0.566546153, -2.03121442, -3.974630891], 
            "desc" => "EDL933_acnB_Z0128  DB_XREF=GI:12512830  SEG=NC_002655:+135988,138585  LEN=2597  DEF=aconitate hydrase B"
        ),
    ]

    # [6] ACt2r::[Acetate reversible transport via proton symport] Acetate reversible transport via proton symport
    # (-1.0) h_e::[H+] + (-1.0) ac_e::[Acetate] <==> (1.0) h_c::[H+] + (1.0) ac_c::[Acetate]
    raw_map["Supp1"]["ACt2r"] = [] # Not find
    raw_map["Supp3"]["ACt2r"] = [] # Not find

    # [7] ADK1::[Adenylate kinase] Adenylate kinase
    # (-1.0) amp_c::[AMP C10H12N5O7P] + (-1.0) atp_c::[ATP C10H12N5O13P3] <==> (2.0) adp_c::[ADP C10H12N5O10P2]
    raw_map["Supp1"]["ADK1"] = ["adk"]
    raw_map["Supp3"]["ADK1"] = [
        Dict{String, Any}(
            "raw" => [7886.2, 9061.0, 9061.2, 10879.3, 9168.4, 5144.6, 7225.0, 3238.9, 7251.7, 3814.0, 895.0, 611.5], 
            "norm" => [0.745125714, 0.945465695, 0.945497538, 1.209309243, 0.962465403, 0.128854325, 0.618793003, -0.538700659, 0.624114658, -0.302899741, -2.394244997, -2.943780181], 
            "desc" => "MG1655_adk_b0474  DB_XREF=PID:g1786680  SEG=NC_000913:+496399,497043  LEN=644  DEF=adenylate kinase activity; pleiotropic effects on glycerol-3-phosphate acyltransferase activity"
        ),
    ]

    # [8] AKGDH::[2-Oxogluterate dehydrogenase] 2-Oxogluterate dehydrogenase
    # (-1.0) nad_c::[Nicotinamide adenine dinucleotide] + (-1.0) akg_c::[2-Oxoglutarate] + (-1.0) coa_c::[Coenzyme A] ==> (1.0) nadh_c::[Nicotinamide adenine dinucleotide - reduced] + (1.0) succoa_c::[Succinyl-CoA] + (1.0) co2_c::[CO2 CO2]
    raw_map["Supp1"]["AKGDH"] = ["sucAB", "lpdA"]
    raw_map["Supp3"]["AKGDH"] = [
        Dict{String, Any}(
            "raw" => [13646.1, 16078.3, 14501.5, 17667.9, 16382.5, 13625.9, 12624.3, 6149.4, 7986.4, 8169.9, 3066.4, 1511.2], 
            "norm" => [0.595442346, 0.832068527, 0.68315579, 0.968084225, 0.859109184, 0.593305177, 0.483157047, -0.554528789, -0.177429111, -0.144656022, -1.558428539, -2.579279835], 
            "desc" => "MG1655_lpdA_b0116  DB_XREF=PID:g1786307  SEG=NC_000913:+127912,129336  LEN=1424  DEF=lipoamide dehydrogenase (NADH); component of 2-oxodehydrogenase and pyruvate complexes; L-protein of glycine cleavage complex"
        ),
    ]

    # [9] AKGt2r::[2 oxoglutarate reversible transport via symport] 2 oxoglutarate reversible transport via symport
    # (-1.0) h_e::[H+] + (-1.0) akg_e::[2-Oxoglutarate] <==> (1.0) h_c::[H+] + (1.0) akg_c::[2-Oxoglutarate]
    raw_map["Supp1"]["AKGt2r"] = [] # Not found
    raw_map["Supp3"]["AKGt2r"] = [] # Not found

    # [10] ALCD2x::[Alcohol dehydrogenase (ethanol)] Alcohol dehydrogenase (ethanol)
    # (-1.0) nad_c::[Nicotinamide adenine dinucleotide] + (-1.0) etoh_c::[Ethanol] <==> (1.0) h_c::[H+] + (1.0) nadh_c::[Nicotinamide adenine dinucleotide - reduced] + (1.0) acald_c::[Acetaldehyde]
    raw_map["Supp1"]["ALCD2x"] = [] # Not Found
    raw_map["Supp3"]["ALCD2x"] = [] # Not Found

    # [11] AMALT1::[Amylomaltase (maltotriose)] Alternate Carbon Metabolism
    # (-1.0) malt_c::[Maltose C12H22O11] ==> (2.0) glc__D_e::[D-Glucose]
    raw_map["Supp1"]["AMALT1"] = [] # Not Found
    raw_map["Supp3"]["AMALT1"] = [
        Dict{String, Any}(
            "raw" => [696.9, 731.6, 781.4, 687.5, 7765.8, 5975.4, 6386.1, 1476.2, 529.0, 406.1, 914.6, 981.8], 
            "norm" => [-0.898532187, -0.828428766, -0.733422586, -0.918124128, 2.579578806, 2.201479545, 2.297379398, 0.184332449, -1.296216119, -1.677648814, -0.506342922, -0.404054675], 
            "desc" => "MG1655_malQ_b3416  DB_XREF=PID:g1789821  SEG=NC_000913:-3545619,3547703  LEN=2084  DEF=4-alpha-glucanotransferase (amylomaltase)"
        ), 
    ]

    # [12] ATPM::[ATP maintenance requirement] ATP maintenance requirement
    # (-1.0) h2o_c::[H2O H2O] + (-1.0) atp_c::[ATP C10H12N5O13P3] ==> (1.0) h_c::[H+] + (1.0) pi_c::[Phosphate] + (1.0) adp_c::[ADP C10H12N5O10P2]
    raw_map["Supp1"]["ATPM"] = [] # Not apply
    raw_map["Supp3"]["ATPM"] = [] # Not apply

    # [13] ATPS4r::[ATP synthase (four protons for one ATP)] ATP synthase (four protons for one ATP)
    # (-4.0) h_e::[H+] + (-1.0) pi_c::[Phosphate] + (-1.0) adp_c::[ADP C10H12N5O10P2] <==> (1.0) h2o_c::[H2O H2O] + (3.0) h_c::[H+] + (1.0) atp_c::[ATP C10H12N5O13P3]
    raw_map["Supp1"]["ATPS4r"] = ["atpABCDEFGHI"]
    raw_map["Supp3"]["ATPS4r"] = [
        Dict{String, Any}(
            "raw" => [12087.3, 12914.6, 12987.4, 13768.2, 13205.6, 9690.6, 10070.9, 4570.8, 9420.4, 7609.0, 3007.2, 1117.5], 
            "norm" => [0.661663269, 0.75717421, 0.765283891, 0.849511209, 0.789321102, 0.342829149, 0.398363867, -0.741310151, 0.302031475, -0.006049983, -1.345336024, -2.773482014], 
            "desc" => "MG1655_atpA_b3734  DB_XREF=PID:g1790172  SEG=NC_000913:-3915944,3917485  LEN=1541  DEF=membrane-bound ATP synthase, F1 sector, alpha-subunit"
        ),
        Dict{String, Any}(
            "raw" => [10510.2, 11655.7, 11213.1, 12381.7, 11434.2, 7936.3, 9853.5, 5822.6, 9795.2, 10527.4, 2942.9, 1054.8], 
            "norm" => [0.500189649, 0.649435176, 0.59358471, 0.736608935, 0.621754957, 0.094937993, 0.407107698, -0.351865057, 0.398546381, 0.502548698, -1.336290049, -2.816559092], 
            "desc" => "MG1655_atpB_b3738  DB_XREF=PID:g1790176  SEG=NC_000913:-3918864,3919679  LEN=815  DEF=membrane-bound ATP synthase, F0 sector, subunit a"
        ),
        Dict{String, Any}(
            "raw" => [4255.3, 4239.2, 3871.9, 5734.1, 5412.9, 4102.9, 4192.1, 2326.0, 4163.3, 5741.9, 2052.2, 908.5], 
            "norm" => [0.265742756, 0.260273943, 0.129523603, 0.696048975, 0.612883646, 0.213125902, 0.244155042, -0.605666993, 0.234209431, 0.698010116, -0.786346752, -1.96195967], 
            "desc" => "MG1655_atpC_b3731  DB_XREF=PID:g1790169  SEG=NC_000913:-3913181,3913600  LEN=419  DEF=membrane-bound ATP synthase, F1 sector, epsilon-subunit"
        ),
        Dict{String, Any}(
            "raw" => [9049.5, 8810.0, 7821.7, 8130.2, 8956.6, 6480.2, 6084.6, 2372.1, 5169.5, 5399.2, 1611.6, 570.1], 
            "norm" => [0.952953004, 0.914256941, 0.742597124, 0.798405764, 0.938066099, 0.471153261, 0.380277345, -0.978720248, 0.14513967, 0.20786058, -1.536391368, -3.035598172], 
            "desc" => "MG1655_atpD_b3732  DB_XREF=PID:g1790170  SEG=NC_000913:-3913621,3915003  LEN=1382  DEF=membrane-bound ATP synthase, F1 sector, beta-subunit"
        ),
        Dict{String, Any}(
            "raw" => [12176.6, 13513.0, 12751.8, 13614.9, 13587.5, 8984.3, 12037.7, 6437.6, 11338.0, 8719.9, 3603.9, 1256.7], 
            "norm" => [0.53801535, 0.688251996, 0.604604903, 0.699090382, 0.696184031, 0.099382003, 0.521463764, -0.381501161, 0.43507017, 0.056287491, -1.218465119, -2.738383809], 
            "desc" => "MG1655_atpE_b3737  DB_XREF=PID:g1790175  SEG=NC_000913:-3918578,3918817  LEN=239  DEF=membrane-bound ATP synthase, F0 sector, subunit c"
        ),
        Dict{String, Any}(
            "raw" => [13639.9, 15432.4, 15390.3, 17036.3, 14945.8, 10935.4, 12808.9, 7634.2, 12064.5, 11461.6, 3565.0, 1180.3], 
            "norm" => [0.528029187, 0.706158563, 0.702217474, 0.84880816, 0.659936242, 0.209202113, 0.437342706, -0.309254992, 0.350964246, 0.277004574, -1.407829898, -3.002578375], 
            "desc" => "MG1655_atpF_b3736  DB_XREF=PID:g1790174  SEG=NC_000913:-3918046,3918516  LEN=470  DEF=membrane-bound ATP synthase, F0 sector, subunit b"
        ),
        Dict{String, Any}(
            "raw" => [10236.4, 8963.2, 9423.6, 9376.8, 10365.9, 6039.8, 7499.6, 2977.1, 7473.5, 5778.6, 1878.4, 584.5],
            "norm" => [0.938861252, 0.747238617, 0.819503031, 0.81232039, 0.956998204, 0.177725506, 0.490038378, -0.84286759, 0.485008775, 0.113944737, -1.507270958, -3.191500342], 
            "desc" => "MG1655_atpG_b3733  DB_XREF=PID:g1790171  SEG=NC_000913:-3915030,3915893  LEN=863  DEF=membrane-bound ATP synthase, F1 sector, gamma-subunit"
        ),
        Dict{String, Any}(
            "raw" => [10769.5, 10578.8, 10359.9, 10599.8, 11439.2, 7385.0, 8407.5, 3590.5, 8224.5, 5064.9, 2155.5, 676.0], 
            "norm" => [0.857976029, 0.832200743, 0.802034835, 0.835061802, 0.945010919, 0.313694584, 0.500773537, -0.726718574, 0.469024637, -0.230369552, -1.462880774, -3.135808185], 
            "desc" => "MG1655_atpH_b3735  DB_XREF=PID:g1790173  SEG=NC_000913:-3917498,3918031  LEN=533  DEF=membrane-bound ATP synthase, F1 sector, delta-subunit"
        ),
        Dict{String, Any}(
            "raw" => [3209.1, 3567.8, 3861.0, 4034.9, 4111.4, 2315.3, 3184.7, 1738.6, 3527.9, 3249.3, 1413.2, 521.1], 
            "norm" => [0.325822097, 0.478688094, 0.592627904, 0.656186268, 0.683283089, -0.14514751, 0.314810826, -0.5584206, 0.462463018, 0.3437823, -0.857380996, -2.29671449], 
            "desc" => "MG1655_atpI_b3739  DB_XREF=PID:g1790177  SEG=NC_000913:-3919688,3920080  LEN=392  DEF=membrane-bound ATP synthase, dispensable protein, affects expression of atpB"
        ),
    ]

    # [14] BIOMASS_Ecoli_core_w_GAM::[Biomass Objective Function with GAM] Biomass Objective Function with GAM
    # (-0.2557) gln__L_c::[L-Glutamine] + (-4.9414) glu__L_c::[L-Glutamate] + (-59.81) h2o_c::[H2O H2O] + (-3.547) nad_c::[Nicotinamide adenine dinucleotide] + (-13.0279) nadph_c::[Nicotinamide adenine dinucleotide phosphate - reduced] + (-1.496) 3pg_c::[3-Phospho-D-glycerate] + (-1.7867) oaa_c::[Oxaloacetate] + (-0.5191) pep_c::[Phosphoenolpyruvate] + (-2.8328) pyr_c::[Pyruvate] + (-0.8977) r5p_c::[Alpha-D-Ribose 5-phosphate] + (-3.7478) accoa_c::[Acetyl-CoA] + (-59.81) atp_c::[ATP C10H12N5O13P3] + (-0.361) e4p_c::[D-Erythrose 4-phosphate] + (-0.0709) f6p_c::[D-Fructose 6-phosphate] + (-0.129) g3p_c::[Glyceraldehyde 3-phosphate] + (-0.205) g6p_c::[D-Glucose 6-phosphate] ==> (59.81) h_c::[H+] + (3.547) nadh_c::[Nicotinamide adenine dinucleotide - reduced] + (13.0279) nadp_c::[Nicotinamide adenine dinucleotide phosphate] + (59.81) pi_c::[Phosphate] + (59.81) adp_c::[ADP C10H12N5O10P2] + (4.1182) akg_c::[2-Oxoglutarate] + (3.7478) coa_c::[Coenzyme A]
    raw_map["Supp1"]["BIOMASS_Ecoli_core_w_GAM"] = []
    raw_map["Supp3"]["BIOMASS_Ecoli_core_w_GAM"] = []

    # [15] CO2t::[CO2 transporter via diffusion] CO2 transporter via diffusion
    # (-1.0) co2_e::[CO2 CO2] <==> (1.0) co2_c::[CO2 CO2]
    raw_map["Supp1"]["CO2t"] = [] # Not apply
    raw_map["Supp3"]["CO2t"] = [] # Not apply

    # [16] CS::[Citrate synthase] Citrate synthase
    # (-1.0) h2o_c::[H2O H2O] + (-1.0) oaa_c::[Oxaloacetate] + (-1.0) accoa_c::[Acetyl-CoA] ==> (1.0) h_c::[H+] + (1.0) cit_c::[Citrate] + (1.0) coa_c::[Coenzyme A]
    raw_map["Supp1"]["CS"] = ["gltA"]
    raw_map["Supp3"]["CS"] = [
        Dict{String, Any}(
            "raw" => [1532.8, 1669.8, 1388.4, 2499.0, 2118.6, 4589.6, 2825.8, 2068.4, 2653.8, 1803.5, 515.8, 483.1], 
            "norm" => [-0.15328134, -0.029775492, -0.296027536, 0.551900095, 0.31366042, 1.428917616, 0.729208554, 0.279064403, 0.638608841, 0.081348616, -1.724567128, -1.819057049], 
            "desc" => "MG1655_gltA_b0720  DB_XREF=PID:g1786939  SEG=NC_000913:-752408,753691  LEN=1283  DEF=citrate synthase"
        ),
        # Dict{String, Any}(
        #     "raw" => [75.3, 79.7, 181.9, 85.5, 48.6, 73.2, 115.9, 228.9, 559.8, 1080.9, 8261.6, 9969.9], 
        #     "norm" => [-1.971277075, -1.889347216, -0.698853302, -1.78800252, -2.602970626, -2.012083292, -1.349118279, -0.367281382, 0.922912642, 1.872162307, 4.806350461, 5.077508284], 
        #     "desc" => "MG1655_prpC_b0333  DB_XREF=PID:g1786527  SEG=NC_000913:+349236,350405  LEN=1169  DEF=putative citrate synthase; propionate metabolism?"
        # ),
    ]

    # [17] CYTBD::[Cytochrome oxidase bd (ubiquinol-8: 2 protons)] Cytochrome oxidase bd (ubiquinol-8: 2 protons)
    # (-2.0) h_c::[H+] + (-0.5) o2_c::[O2 O2] + (-1.0) q8h2_c::[Ubiquinol-8] ==> (1.0) h2o_c::[H2O H2O] + (2.0) h_e::[H+] + (1.0) q8_c::[Ubiquinone-8]
    raw_map["Supp1"]["CYTBD"] = ["cydABCD", "appBC", "cycBC"]
    raw_map["Supp3"]["CYTBD"] = [
        Dict{String, Any}(
            "raw" => [3344.4, 3200.5, 2862.1, 3087.0, 6081.9, 4348.9, 2944.1, 3392.5, 3485.3, 4036.5, 3868.4, 4079.0], 
            "norm" => [-0.12517031, -0.18862041, -0.349843639, -0.240712236, 0.737604377, 0.253732817, -0.309091043, -0.104568903, -0.065634876, 0.146187173, 0.084819261, 0.161297789], 
            "desc" => "MG1655_cydA_b0733  DB_XREF=PID:g1786953  SEG=NC_000913:+770678,772249  LEN=1571  DEF=cytochrome d terminal oxidase, polypeptide subunit I"
        ),
        Dict{String, Any}(
            "raw" => [2490.4, 2257.9, 2384.8, 2754.1, 5153.9, 3256.7, 2489.9, 3262.3, 2893.4, 4507.6, 3307.0, 4110.9], 
            "norm" => [-0.330952813, -0.472348703, -0.393462014, -0.185749351, 0.718334251, 0.056080533, -0.331242493, 0.058559163, -0.114564511, 0.525029203, 0.078192751, 0.392123983], 
            "desc" => "MG1655_cydB_b0734  DB_XREF=PID:g1786954  SEG=NC_000913:+772265,773404  LEN=1139  DEF=cytochrome d terminal oxidase polypeptide subunit II"
        ),
        Dict{String, Any}(
            "raw" => [686.9, 558.9, 641.1, 471.8, 688.1, 351.6, 438.8, 83.4, 351.3, 215.2, 91.3, 119.8], 
            "norm" => [1.138654456, 0.841144546, 1.03910378, 0.59672979, 1.141172615, 0.172489442, 0.492117897, -1.90332634, 0.171257948, -0.535767551, -1.772758863, -1.38081772], 
            "desc" => "MG1655_cydC_b0886  DB_XREF=PID:g1787112  SEG=NC_000913:-926697,928418  LEN=1721  DEF=ATP-binding component of cytochrome-related transport"
        ),
        Dict{String, Any}(
            "raw" => [544.1, 481.8, 604.1, 450.5, 641.0, 310.0, 350.0, 107.9, 391.2, 145.8, 137.5, 108.5], 
            "norm" => [0.87743126, 0.701993825, 1.028346818, 0.605086537, 1.113883788, 0.065827647, 0.240914353, -1.456745704, 0.401465802, -1.022449849, -1.10700895, -1.448745526], 
            "desc" => "MG1655_cydD_b0887  DB_XREF=PID:g1787113  SEG=NC_000913:-928419,930185  LEN=1766  DEF=ATP-binding component of cytochrome-related transport, Zn sensitive"
        ),

        Dict{String, Any}(
            "raw" => [270.0, 230.5, 253.5, 206.5, 203.6, 276.4, 256.0, 525.3, 407.6, 205.2, 558.0, 807.1],
            "norm" => [-0.219109893, -0.447302549, -0.310083553, -0.605927518, -0.626331738, -0.185311684, -0.29592549, 0.741072285, 0.375084752, -0.615038569, 0.828195822, 1.360678135],
            "desc" => "MG1655_appB_b0979  DB_XREF=PID:g1787213  SEG=NC_000913:+1038519,1039655  LEN=1136  DEF=probable third cytochrome oxidase, subunit II"
        ),
        Dict{String, Any}(
            "raw" => [436.2, 350.3, 312.7, 224.7, 232.7, 281.6, 312.3, 621.9, 489.4, 250.1, 757.6, 900.6],
            "norm" => [0.171328453, -0.145070329, -0.308882097, -0.785661192, -0.735190106, -0.460013983, -0.310728747, 0.6830213, 0.337352785, -0.631156259, 0.967775014, 1.217225161],
            "desc" => "MG1655_appC_b0978  DB_XREF=PID:g1787212  SEG=NC_000913:+1036963,1038507  LEN=1544  DEF=probable third cytochrome oxidase, subunit I"
        ),

        Dict{String, Any}(
            "raw" => [11405.0, 13045.9, 12950.4, 15057.8, 13476.6, 6472.3, 9090.9, 3286.8, 9087.9, 7695.6, 2157.3, 1218.7],
            "norm" => [0.703694248, 0.897624276, 0.887024461, 1.104538804, 0.944484367, -0.113621813, 0.376522835, -1.091216622, 0.376046666, 0.136133519, -1.698673478, -2.522557263],
            "desc" => "MG1655_cyoA_b0432  DB_XREF=PID:g1786635  SEG=NC_000913:-449887,450834  LEN=947  DEF=cytochrome o ubiquinol oxidase subunit II"
        ),
        Dict{String, Any}(
            "raw" => [11036.7, 12868.3, 12488.7, 13469.9, 11319.8, 5924.0, 6449.2, 1962.0, 6394.7, 4177.8, 1106.2, 433.6],
            "norm" => [1.134847992, 1.3563606, 1.313162433, 1.422278265, 1.171387593, 0.23718267, 0.35973124, -1.357063929, 0.347487708, -0.266645541, -2.183776723, -3.534952309],
            "desc" => "MG1655_cyoB_b0431  DB_XREF=PID:g1786634  SEG=NC_000913:-447874,449865  LEN=1991  DEF=cytochrome o ubiquinol oxidase subunit I"
        ),
        Dict{String, Any}(
            "raw" => [9258.5, 9642.8, 9960.4, 11057.3, 10072.2, 4826.7, 5838.4, 1752.9, 5429.8, 6814.7, 1373.3, 531.0],
            "norm" => [0.962778297, 1.021451946, 1.068203501, 1.218927063, 1.08430675, 0.02303698, 0.297572876, -1.438256485, 0.192898879, 0.520649967, -1.790353361, -3.161216414],
            "desc" => "MG1655_cyoC_b0430  DB_XREF=PID:g1786633  SEG=NC_000913:-447270,447884  LEN=614  DEF=cytochrome o ubiquinol oxidase subunit III"
        ),
        Dict{String, Any}(
            "raw" => [8274.3, 9010.5, 8514.5, 9656.1, 9265.9, 4590.2, 5381.9, 1892.2, 5258.0, 8275.6, 2183.0, 918.2],
            "norm" => [0.733066793, 0.856036691, 0.774351339, 0.955870143, 0.896360639, -0.117013459, 0.112545111, -1.395505888, 0.078943669, 0.733293441, -1.189258343, -2.438690136],
            "desc" => "MG1655_cyoD_b0429  DB_XREF=PID:g1786632  SEG=NC_000913:-446941,447270  LEN=329  DEF=cytochrome o ubiquinol oxidase subunit IV"
        ),
        Dict{String, Any}(
            "raw" => [4913.1, 5335.2, 5416.2, 6667.9, 6168.8, 3158.1, 4114.6, 1705.5, 3903.4, 6035.7, 1876.2, 971.4],
            "norm" => [0.429476392, 0.548385143, 0.570123802, 0.870075257, 0.757832662, -0.208100357, 0.173594975, -1.096962456, 0.097574099, 0.726363888, -0.959343586, -1.90901982],
            "desc" => "MG1655_cyoE_b0428  DB_XREF=PID:g1786631  SEG=NC_000913:-446039,446929  LEN=890  DEF=protoheme IX farnesyltransferase (haeme O biosynthesis)"
        ),
    ]

    # [18] D_LACt2::[D lactate transport via proton symport] D lactate transport via proton symport
    # (-1.0) h_e::[H+] + (-1.0) lac__D_e::[D-Lactate] <==> (1.0) h_c::[H+] + (1.0) lac__D_c::[D-Lactate]
    raw_map["Supp1"]["D_LACt2"] = [] # Not found
    raw_map["Supp3"]["D_LACt2"] = [] # Not found

    # [19] ENO::[Enolase] Enolase
    # (-1.0) 2pg_c::[D-Glycerate 2-phosphate] <==> (1.0) h2o_c::[H2O H2O] + (1.0) pep_c::[Phosphoenolpyruvate]
    raw_map["Supp1"]["ENO"] = []
    raw_map["Supp3"]["ENO"] = [
        Dict{String, Any}(
            "raw" => [10503.6, 13024.1, 12521.9, 13701.7, 12662.8, 12241.3, 16149.4, 13807.2, 12180.4, 10351.9, 18905.6, 20348.1], 
            "norm" => [-0.370640731, -0.06034093, -0.117071127, 0.01283029, -0.100928163, -0.149767835, 0.249955953, 0.023896169, -0.156963101, -0.391629026, 0.477289024, 0.583369477], 
            "desc" => "MG1655_eno_b2779  DB_XREF=PID:g1789141  SEG=NC_000913:-2904665,2905963  LEN=1298  DEF=enolase"
        ),
    ]

    # [20] ETOHt2r::[Ethanol reversible transport via proton symport] Ethanol reversible transport via proton symport
    # (-1.0) h_e::[H+] + (-1.0) etoh_e::[Ethanol] <==> (1.0) h_c::[H+] + (1.0) etoh_c::[Ethanol]
    raw_map["Supp1"]["ETOHt2r"] = [] # Not found
    raw_map["Supp3"]["ETOHt2r"] = [] # Not found

    # [21] EX_ac_e::[Acetate exchange] Acetate exchange
    # (-1.0) ac_e::[Acetate] <==> 
    raw_map["Supp1"]["EX_ac_e"] = [] # Not apply
    raw_map["Supp3"]["EX_ac_e"] = [] # Not apply

    # [22] EX_acald_e::[Acetaldehyde exchange] Acetaldehyde exchange
    # (-1.0) acald_e::[Acetaldehyde] ==> 
    raw_map["Supp1"]["EX_acald_e"] = [] # Not apply
    raw_map["Supp3"]["EX_acald_e"] = [] # Not apply

    # [23] EX_akg_e::[2-Oxoglutarate exchange] 2-Oxoglutarate exchange
    # (-1.0) akg_e::[2-Oxoglutarate] ==> 
    raw_map["Supp1"]["EX_akg_e"] = [] # Not apply
    raw_map["Supp3"]["EX_akg_e"] = [] # Not apply

    # [24] EX_co2_e::[CO2 exchange] CO2 exchange
    # (-1.0) co2_e::[CO2 CO2] <==> 
    raw_map["Supp1"]["EX_co2_e"] = [] # Not apply
    raw_map["Supp3"]["EX_co2_e"] = [] # Not apply

    # [25] EX_etoh_e::[Ethanol exchange] Ethanol exchange
    # (-1.0) etoh_e::[Ethanol] ==> 
    raw_map["Supp1"]["EX_etoh_e"] = []
    raw_map["Supp3"]["EX_etoh_e"] = []

    # [26] EX_for_e::[Formate exchange] Formate exchange
    # (-1.0) for_e::[Formate] ==> 
    raw_map["Supp1"]["EX_for_e"] = [] # Not apply
    raw_map["Supp3"]["EX_for_e"] = [] # Not apply

    # [27] EX_fru_e::[D-Fructose exchange] D-Fructose exchange
    # (-1.0) fru_e::[D-Fructose] ==> 
    raw_map["Supp1"]["EX_fru_e"] = [] # Not apply
    raw_map["Supp3"]["EX_fru_e"] = [] # Not apply

    # [28] EX_fum_e::[Fumarate exchange] Fumarate exchange
    # (-1.0) fum_e::[Fumarate] ==> 
    raw_map["Supp1"]["EX_fum_e"] = [] # Not apply
    raw_map["Supp3"]["EX_fum_e"] = [] # Not apply

    # [29] EX_gal_e::[D-Galactose exchange] Extracellular exchange
    # (-1.0) gal_e::[D-Galactose] <==> 
    raw_map["Supp1"]["EX_gal_e"] = [] # Not apply
    raw_map["Supp3"]["EX_gal_e"] = [] # Not apply

    # [30] EX_glc__D_e::[D-Glucose exchange] D-Glucose exchange
    # (-1.0) glc__D_e::[D-Glucose] <==> 
    raw_map["Supp1"]["EX_glc__D_e"] = [] # Not apply
    raw_map["Supp3"]["EX_glc__D_e"] = [] # Not apply

    # [31] EX_gln__L_e::[L-Glutamine exchange] L-Glutamine exchange
    # (-1.0) gln__L_e::[L-Glutamine] ==> 
    raw_map["Supp1"]["EX_gln__L_e"] = [] # Not apply
    raw_map["Supp3"]["EX_gln__L_e"] = [] # Not apply

    # [32] EX_glu__L_e::[L-Glutamate exchange] L-Glutamate exchange
    # (-1.0) glu__L_e::[L-Glutamate] ==> 
    raw_map["Supp1"]["EX_glu__L_e"] = [] # Not apply
    raw_map["Supp3"]["EX_glu__L_e"] = [] # Not apply

    # [33] EX_glyc_e::[Glycerol exchange] Extracellular exchange
    # (-1.0) glyc_e::[Glycerol] <==> 
    raw_map["Supp1"]["EX_glyc_e"] = [] # Not apply
    raw_map["Supp3"]["EX_glyc_e"] = [] # Not apply

    # [34] EX_h2o_e::[H2O exchange] H2O exchange
    # (-1.0) h2o_e::[H2O H2O] <==> 
    raw_map["Supp1"]["EX_h2o_e"] = [] # Not apply
    raw_map["Supp3"]["EX_h2o_e"] = [] # Not apply

    # [35] EX_h_e::[H+ exchange] H+ exchange
    # (-1.0) h_e::[H+] <==> 
    raw_map["Supp1"]["EX_h_e"] = [] # Not apply
    raw_map["Supp3"]["EX_h_e"] = [] # Not apply

    # [36] EX_lac__D_e::[D-lactate exchange] D-lactate exchange
    # (-1.0) lac__D_e::[D-Lactate] <==> 
    raw_map["Supp1"]["EX_lac__D_e"] = [] # Not apply
    raw_map["Supp3"]["EX_lac__D_e"] = [] # Not apply

    # [37] EX_mal__L_e::[L-Malate exchange] L-Malate exchange
    # (-1.0) mal__L_e::[L-Malate] ==> 
    raw_map["Supp1"]["EX_mal__L_e"] = [] # Not apply
    raw_map["Supp3"]["EX_mal__L_e"] = [] # Not apply

    # [38] EX_malt_e::[Maltose exchange] Extracellular exchange
    # (-1.0) malt_e::[Maltose C12H22O11] <==> 
    raw_map["Supp1"]["EX_malt_e"] = [] # Not apply
    raw_map["Supp3"]["EX_malt_e"] = [] # Not apply

    # [39] EX_nh4_e::[Ammonia exchange] Ammonia exchange
    # (-1.0) nh4_e::[Ammonium] <==> 
    raw_map["Supp1"]["EX_nh4_e"] = [] # Not apply
    raw_map["Supp3"]["EX_nh4_e"] = [] # Not apply

    # [40] EX_o2_e::[O2 exchange] O2 exchange
    # (-1.0) o2_e::[O2 O2] <==> 
    raw_map["Supp1"]["EX_o2_e"] = [] # Not apply
    raw_map["Supp3"]["EX_o2_e"] = [] # Not apply

    # [41] EX_pi_e::[Phosphate exchange] Phosphate exchange
    # (-1.0) pi_e::[Phosphate] <==> 
    raw_map["Supp1"]["EX_pi_e"] = [] # Not apply
    raw_map["Supp3"]["EX_pi_e"] = [] # Not apply

    # [42] EX_pyr_e::[Pyruvate exchange] Pyruvate exchange
    # (-1.0) pyr_e::[Pyruvate] ==> 
    raw_map["Supp1"]["EX_pyr_e"] = [] # Not apply
    raw_map["Supp3"]["EX_pyr_e"] = [] # Not apply

    # [43] EX_succ_e::[Succinate exchange] Succinate exchange
    # (-1.0) succ_e::[Succinate] ==> 
    raw_map["Supp1"]["EX_succ_e"] = [] # Not apply
    raw_map["Supp3"]["EX_succ_e"] = [] # Not apply

    # [44] FBA::[Fructose-bisphosphate aldolase] Fructose-bisphosphate aldolase
    # (-1.0) fdp_c::[D-Fructose 1,6-bisphosphate] <==> (1.0) dhap_c::[Dihydroxyacetone phosphate] + (1.0) g3p_c::[Glyceraldehyde 3-phosphate]
    raw_map["Supp1"]["FBA"] = ["fbaA"]
    raw_map["Supp3"]["FBA"] = [
        Dict{String, Any}(
            "raw" => [11465.6, 12410.1, 13292.5, 13914.2, 11092.9, 10992.3, 11961.8, 9545.6, 9329.4, 5453.1, 7897.6, 7774.3], 
            "norm" => [0.181083122, 0.295286009, 0.394383734, 0.460329231, 0.133407844, 0.120264551, 0.242205769, -0.083320944, -0.116372527, -0.891080215, -0.356742527, -0.379444047], 
            "desc" => "MG1655_b2925  DB_XREF=PID:g1789293  SEG=NC_000913:-3068185,3069264  LEN=1079  DEF=fructose-bisphosphate aldolase, class II"
        ),
    ]

    # [45] FBP::[Fructose-bisphosphatase] Fructose-bisphosphatase
    # (-1.0) h2o_c::[H2O H2O] + (-1.0) fdp_c::[D-Fructose 1,6-bisphosphate] ==> (1.0) pi_c::[Phosphate] + (1.0) f6p_c::[D-Fructose 6-phosphate]
    raw_map["Supp1"]["FBP"] = ["fbp"]
    raw_map["Supp3"]["FBP"] = [
        Dict{String, Any}(
            "raw" => [1796.1, 1893.9, 2315.9, 2053.4, 1816.1, 1391.1, 1127.7, 1111.1, 1522.7, 1075.3, 671.1, 609.8], 
            "norm" => [0.419814683, 0.496307164, 0.786519967, 0.612961697, 0.435790651, 0.05117314, -0.251669672, -0.273064326, 0.18157874, -0.320313777, -1.000453331, -1.138644937], 
            "desc" => "MG1655_fbp_b4232  DB_XREF=PID:g1790679  SEG=NC_000913:-4452185,4453183  LEN=998  DEF=fructose-bisphosphatase"
        ),
    ]

    # [46] FORt::[Formate transport via diffusion] Formate transport via diffusion
    # (-1.0) for_e::[Formate] <== (1.0) for_c::[Formate]
    raw_map["Supp1"]["FORt"] = [] # Not apply
    raw_map["Supp3"]["FORt"] = [] # Not apply

    # [47] FORt2::[Formate transport in via proton symport] Formate transport in via proton symport
    # (-1.0) h_e::[H+] + (-1.0) for_e::[Formate] ==> (1.0) h_c::[H+] + (1.0) for_c::[Formate]
    raw_map["Supp1"]["FORt2"] = [] # Not found
    raw_map["Supp3"]["FORt2"] = [] # Not found

    # [48] FRD7::[Fumarate reductase] Fumarate reductase
    # (-1.0) q8h2_c::[Ubiquinol-8] + (-1.0) fum_c::[Fumarate] ==> (1.0) q8_c::[Ubiquinone-8] + (1.0) succ_c::[Succinate]
    raw_map["Supp1"]["FRD7"] = ["frdABCD"]
    raw_map["Supp3"]["FRD7"] = [
        Dict{String, Any}(
            "raw" => [520.9, 472.1, 574.3, 651.8, 921.9, 828.8, 603.0, 1568.7, 507.9, 1108.6, 1698.4, 1905.0], 
            "norm" => [-0.687215031, -0.829128986, -0.546416907, -0.363792117, 0.136388799, -0.017197466, -0.476063466, 0.903276102, -0.723676995, 0.402445539, 1.017882903, 1.183497624], 
            "desc" => "MG1655_frdA_b4154  DB_XREF=PID:g1790597  SEG=NC_000913:-4378088,4379896  LEN=1808  DEF=fumarate reductase, anaerobic, flavoprotein subunit"
        ),
        Dict{String, Any}(
            "raw" => [437.0, 416.1, 496.8, 557.0, 814.9, 807.5, 681.5, 1736.9, 562.9, 1365.6, 2384.6, 2469.7], 
            "norm" => [-0.986402453, -1.057105444, -0.801370559, -0.636358405, -0.087412702, -0.100573473, -0.345322076, 1.004407057, -0.621157085, 0.657427326, 1.461639647, 1.512228167], 
            "desc" => "MG1655_frdB_b4153  DB_XREF=PID:g1790596  SEG=NC_000913:-4377361,4378095  LEN=734  DEF=fumarate reductase, anaerobic, iron-sulfur protein subunit"
        ),
        Dict{String, Any}(
            "raw" => [173.0, 227.3, 285.8, 305.6, 339.0, 484.1, 349.4, 752.8, 274.1, 687.6, 847.4, 986.5], 
            "norm" => [-1.249010532, -0.855184886, -0.524776653, -0.428138026, -0.278497296, 0.235522524, -0.234902962, 0.872484058, -0.585080242, 0.741786975, 1.043260559, 1.262536481], 
            "desc" => "MG1655_frdC_b4152  DB_XREF=PID:g1790595  SEG=NC_000913:-4376955,4377350  LEN=395  DEF=fumarate reductase, anaerobic, membrane anchor polypeptide"
        ),
        Dict{String, Any}(
            "raw" => [333.5, 419.0, 440.4, 451.2, 673.9, 703.2, 674.8, 1496.7, 432.8, 1084.0, 1512.4, 2080.9], 
            "norm" => [-1.111626423, -0.782362941, -0.710498716, -0.675546117, -0.096778659, -0.035378114, -0.094853211, 1.054399986, -0.735612685, 0.588979667, 1.069454665, 1.529822547], 
            "desc" => "MG1655_frdD_b4151  DB_XREF=PID:g1790594  SEG=NC_000913:-4376585,4376944  LEN=359  DEF=fumarate reductase, anaerobic, membrane anchor polypeptide"
        ),
    ]

    # [49] FRUpts2::[Fructose transport via PEP:Pyr PTS (f6p generating)] Fructose transport via PEP:Pyr PTS (f6p generating)
    # (-1.0) pep_c::[Phosphoenolpyruvate] + (-1.0) fru_e::[D-Fructose] ==> (1.0) pyr_c::[Pyruvate] + (1.0) f6p_c::[D-Fructose 6-phosphate]
    raw_map["Supp1"]["FRUpts2"] = [] # Not found
    raw_map["Supp3"]["FRUpts2"] = [] # Not found

    # [50] FUM::[Fumarase] Fumarase
    # (-1.0) h2o_c::[H2O H2O] + (-1.0) fum_c::[Fumarate] <==> (1.0) mal__L_c::[L-Malate]
    raw_map["Supp1"]["FUM"] = ["fumA", "fumB", "fumC"]
    raw_map["Supp3"]["FUM"] = [
        Dict{String, Any}(
            "raw" => [4784.2, 4693.5, 4557.1, 5988.7, 6264.7, 6557.5, 4276.1, 4075.5, 2429.0, 8003.1, 3783.1, 2213.0], 
            "norm" => [0.084156158, 0.056542615, 0.013994485, 0.408121319, 0.47312388, 0.539024359, -0.077825952, -0.147144478, -0.893759054, 0.826437392, -0.254552632, -1.028118093], 
            "desc" => "MG1655_fumA_b1612  DB_XREF=PID:g1787897  SEG=NC_000913:-1684755,1686401  LEN=1646  DEF=fumarase A = fumarate hydratase Class I; aerobic isozyme"
        ),
        Dict{String, Any}(
            "raw" => [93.8, 114.1, 290.4, 104.8, 97.4, 112.2, 86.2, 137.8, 65.3, 55.5, 41.6, 209.1], 
            "norm" => [-0.126890495, 0.155748469, 1.50349113, 0.033088394, -0.072556645, 0.131522353, -0.248790548, 0.428025565, -0.649395426, -0.883990646, -1.299894889, 1.029642739], 
            "desc" => "CFT073_fumB_c5127  DB_XREF=26250935  SEG=NC_004431:-4902575,4904221  LEN=1646  DEF=Fumarate hydratase class I, anaerobic"
        ),
        Dict{String, Any}(
            "raw" => [868.7, 730.7, 643.6, 667.1, 763.6, 1020.5, 536.9, 481.3, 308.0, 792.3, 609.3, 517.3], 
            "norm" => [0.454279129, 0.204700299, 0.021585417, 0.073324133, 0.268238194, 0.686625368, -0.239925504, -0.397642485, -1.041648558, 0.321467893, -0.057426168, -0.293577717], 
            "desc" => "MG1655_fumC_b1611  DB_XREF=PID:g1787896  SEG=NC_000913:-1683209,1684612  LEN=1403  DEF=fumarase C= fumarate hydratase Class II; isozyme"
        ),
    ]

    # [51] FUMt2_2::[Fumarate transport via proton symport (2 H)] Fumarate transport via proton symport (2 H)
    # (-2.0) h_e::[H+] + (-1.0) fum_e::[Fumarate] ==> (2.0) h_c::[H+] + (1.0) fum_c::[Fumarate]
    raw_map["Supp1"]["FUMt2_2"] = []
    raw_map["Supp3"]["FUMt2_2"] = []

    # [52] G3PD2::[Glycerol-3-phosphate dehydrogenase (NADP)] Alternate Carbon Metabolism
    # (-1.0) nadp_c::[Nicotinamide adenine dinucleotide phosphate] + (-1.0) glyc3p_c::[Glycerol 3-phosphate] <==> (1.0) h_c::[H+] + (1.0) nadph_c::[Nicotinamide adenine dinucleotide phosphate - reduced] + (1.0) dhap_c::[Dihydroxyacetone phosphate]
    raw_map["Supp1"]["G3PD2"] = ["gpsA"]
    raw_map["Supp3"]["G3PD2"] = [
        Dict{String, Any}(
            "raw" => [3437.0, 3821.3, 3929.4, 3343.3, 3155.0, 3106.6, 2978.2, 3373.0, 3520.9, 2398.3, 2887.6, 3271.1], 
            "norm" => [0.083662531, 0.236576203, 0.276821716, 0.043785495, -0.039847316, -0.062150826, -0.12304668, 0.056544998, 0.118456932, -0.435475186, -0.167616412, 0.012288544], 
            "desc" => "MG1655_gpsA_b3608  DB_XREF=PID:g1790037  SEG=NC_000913:-3780269,3781288  LEN=1019  DEF=glycerol-3-phosphate dehydrogenase (NAD+)"
        ),
    ]

    # [53] G6PDH2r::[Glucose 6-phosphate dehydrogenase] Glucose 6-phosphate dehydrogenase
    # (-1.0) nadp_c::[Nicotinamide adenine dinucleotide phosphate] + (-1.0) g6p_c::[D-Glucose 6-phosphate] <==> (1.0) h_c::[H+] + (1.0) nadph_c::[Nicotinamide adenine dinucleotide phosphate - reduced] + (1.0) 6pgl_c::[6-phospho-D-glucono-1,5-lactone]
    raw_map["Supp1"]["G6PDH2r"] = ["zwf"]
    raw_map["Supp3"]["G6PDH2r"] = [
        Dict{String, Any}(
            "raw" => [5233.0, 3937.0, 3881.3, 3660.1, 3493.1, 3331.8, 3103.2, 2624.5, 2955.0, 2293.5, 1666.4, 1313.1], 
            "norm" => [0.83345596, 0.422914414, 0.40235765, 0.317700767, 0.250325646, 0.182119503, 0.079574382, -0.162139701, 0.008975832, -0.356631391, -0.817447554, -1.161205509], 
            "desc" => "MG1655_zwf_b1852  DB_XREF=PID:g1788158  SEG=NC_000913:-1932863,1934338  LEN=1475  DEF=glucose-6-phosphate dehydrogenase"
        ),
    ]

    # [54] GALKr::[Galactokinase] Alternate Carbon Metabolism
    # (-1.0) atp_c::[ATP C10H12N5O13P3] + (-1.0) gal_c::[D-Galactose] <==> (1.0) h_c::[H+] + (1.0) adp_c::[ADP C10H12N5O10P2] + (1.0) gal1p_c::[Alpha-D-Galactose 1-phosphate]
    raw_map["Supp1"]["GALKr"] = ["galK"]
    raw_map["Supp3"]["GALKr"] = [
        Dict{String, Any}(
            "raw" => [1272.0, 1335.1, 2775.8, 5023.4, 6883.3, 6595.7, 3769.6, 1332.6, 646.7, 699.2, 998.5, 813.3], 
            "norm" => [-0.565695923, -0.495846789, 0.56010903, 1.415869564, 1.870305795, 1.808731187, 1.00161685, -0.498550795, -1.541626078, -1.429017504, -0.914960261, -1.210935075], 
            "desc" => "MG1655_galK_b0757  DB_XREF=PID:g1786972  SEG=NC_000913:-788054,789202  LEN=1148  DEF=galactokinase"
        ),
    ]

    # [55] GALabcpp::[D-galactose transport via ABC system (periplasm)] Transport, Inner Membrane
    # (-1.0) h2o_e::[H2O H2O] + (-1.0) atp_c::[ATP C10H12N5O13P3] + (-1.0) gal_e::[D-Galactose] ==> (1.0) h_c::[H+] + (1.0) pi_c::[Phosphate] + (1.0) adp_c::[ADP C10H12N5O10P2] + (1.0) gal_c::[D-Galactose]
    raw_map["Supp1"]["GALabcpp"] = []
    raw_map["Supp3"]["GALabcpp"] = [
        Dict{String, Any}(
            "raw" => [193.9, 235.3, 442.4, 666.9, 549.3, 715.6, 375.2, 104.6, 83.8, 64.7, 114.4, 157.1], 
            "norm" => [-0.236733556, 0.042454961, 0.953305026, 1.545424089, 1.265547932, 1.647107028, 0.715613468, -1.127163508, -1.447024211, -1.820208742, -0.997959308, -0.540363179], 
            "desc" => "MG1655_galP_b2943  DB_XREF=PID:g1789312  SEG=NC_000913:+3086303,3087697  LEN=1394  DEF=galactose-proton symport of transport system"
        ),
        Dict{String, Any}(
            "raw" => [12.5, 26.8, 16.4, 26.8, 29.7, 22.6, 21.9, 24.7, 14.3, 15.8, 33.0, 40.8], 
            "norm" => [-0.843846272, 0.256458633, -0.452078553, 0.256458633, 0.404688564, 0.010548405, -0.034843498, 0.138736674, -0.64975922, -0.505849809, 0.556691657, 0.862794785], 
            "desc" => "EDL933_galP_Z4288  DB_XREF=GI:12517486  SEG=NC_002655:+3897091,3898485  LEN=1394  DEF=galactose-proton symport of transport system"
        ),
    ]

    # [56] GAPD::[Glyceraldehyde-3-phosphate dehydrogenase] Glyceraldehyde-3-phosphate dehydrogenase
    # (-1.0) nad_c::[Nicotinamide adenine dinucleotide] + (-1.0) pi_c::[Phosphate] + (-1.0) g3p_c::[Glyceraldehyde 3-phosphate] <==> (1.0) h_c::[H+] + (1.0) nadh_c::[Nicotinamide adenine dinucleotide - reduced] + (1.0) 13dpg_c::[3-Phospho-D-glyceroyl phosphate]
    raw_map["Supp1"]["GAPD"] = ["gapA", "gapC1C2"]
    raw_map["Supp3"]["GAPD"] = [
        Dict{String, Any}(
            "raw" => [16793.4, 20261.8, 19891.5, 24910.6, 16118.8, 17220.1, 26433.1, 23090.0, 20280.1, 19516.5, 18953.0, 19002.2], 
            "norm" => [-0.251386652, 0.019481344, -0.007128928, 0.317478771, -0.310536657, -0.21518748, 0.40306463, 0.207987173, 0.020783766, -0.034586651, -0.076854775, -0.073114542], 
            "desc" => "MG1655_gapA_b1779  DB_XREF=PID:g1788079  SEG=NC_000913:+1860795,1861790  LEN=995  DEF=glyceraldehyde-3-phosphate dehydrogenase A",
        ),
        Dict{String, Any}(
            "raw" => [1282.9, 1972.5, 2188.6, 1795.8, 1729.7, 1739.9, 2079.3, 1868.0, 2613.1, 1198.8, 1166.1, 1414.0], 
            "norm" => [-0.409759162, 0.210857419, 0.360840423, 0.075458804, 0.021353957, 0.029836509, 0.286930043, 0.132326574, 0.616594455, -0.507576892, -0.547476368, -0.269385761], 
            "desc" => "MG1655_gapC_1_b1417  DB_XREF=PID:g1787686  SEG=NC_000913:-1487985,1488389  LEN=404  DEF=glyceraldehyde 3-phosphate dehydrogenase C, interrupted",
        ),
        Dict{String, Any}(
            "raw" => [946.1, 1458.9, 1538.4, 1133.6, 1068.8, 921.7, 1017.1, 648.8, 1312.2, 608.2, 679.5, 854.2], 
            "norm" => [-0.041248152, 0.58356826, 0.660117931, 0.219598926, 0.134679176, -0.078943581, 0.063148793, -0.585467012, 0.430674889, -0.678695016, -0.518767281, -0.188666934], 
            "desc" => "MG1655_gapC_2_b1416  DB_XREF=PID:g1787685  SEG=NC_000913:-1487737,1487988  LEN=251  DEF=glyceraldehyde-3-phosphate dehydrogenase (second fragment)",
        ),
    ]

    # [57] GLCpts::[D-glucose transport via PEP:Pyr PTS] D-glucose transport via PEP:Pyr PTS
    # (-1.0) glc__D_e::[D-Glucose] + (-1.0) pep_c::[Phosphoenolpyruvate] ==> (1.0) pyr_c::[Pyruvate] + (1.0) g6p_c::[D-Glucose 6-phosphate]
    raw_map["Supp1"]["GLCpts"] = []
    raw_map["Supp3"]["GLCpts"] = [
        Dict{String, Any}(
            "raw" => [8898.2, 9431.8, 8300.3, 10712.6, 3524.0, 3969.2, 3978.1, 3602.6, 2944.5, 2780.2, 3600.0, 3062.1], 
            "norm" => [0.90281713, 0.986836731, 0.802467085, 1.170540371, -0.433482472, -0.261848137, -0.258616853, -0.401657919, -0.692673722, -0.775507726, -0.402699489, -0.636174998], 
            "desc" => "MG1655_ptsG_b1101  DB_XREF=PID:g1787343  SEG=NC_000913:+1157092,1158525  LEN=1433  DEF=PTS system, glucose-specific IIBC component",
        ),
    ]

    # [58] GLNS::[Glutamine synthetase] Glutamine synthetase
    # (-1.0) glu__L_c::[L-Glutamate] + (-1.0) nh4_c::[Ammonium] + (-1.0) atp_c::[ATP C10H12N5O13P3] ==> (1.0) gln__L_c::[L-Glutamine] + (1.0) h_c::[H+] + (1.0) pi_c::[Phosphate] + (1.0) adp_c::[ADP C10H12N5O10P2]
    raw_map["Supp1"]["GLNS"] = ["glnA"]
    raw_map["Supp3"]["GLNS"] = [
        Dict{String, Any}(
            "raw" => [7847.6, 8828.6, 7786.1, 8010.7, 7563.1, 4192.7, 4954.9, 718.9, 1747.0, 1013.8, 307.4, 284.0], 
            "norm" => [1.621947657, 1.791880829, 1.610597023, 1.651624465, 1.568673844, 0.717575754, 0.95855209, -1.826440842, -0.545424242, -1.330530782, -3.052114781, -3.166341016], 
            "desc" => "CFT073_glnA_c4819  DB_XREF=26250634  SEG=NC_004431:-4582008,4583417  LEN=1409  DEF=Glutamine synthetase",
        ),
    ]

    # [59] GLNabc::[L-glutamine transport via ABC system] L-glutamine transport via ABC system
    # (-1.0) gln__L_e::[L-Glutamine] + (-1.0) h2o_c::[H2O H2O] + (-1.0) atp_c::[ATP C10H12N5O13P3] ==> (1.0) gln__L_c::[L-Glutamine] + (1.0) h_c::[H+] + (1.0) pi_c::[Phosphate] + (1.0) adp_c::[ADP C10H12N5O10P2]
    raw_map["Supp1"]["GLNabc"] = [] # Not found
    raw_map["Supp3"]["GLNabc"] = [] # Not found

    # [60] GLUDy::[Glutamate dehydrogenase (NADP)] Glutamate dehydrogenase (NADP)
    # (-1.0) glu__L_c::[L-Glutamate] + (-1.0) h2o_c::[H2O H2O] + (-1.0) nadp_c::[Nicotinamide adenine dinucleotide phosphate] <==> (1.0) h_c::[H+] + (1.0) nadph_c::[Nicotinamide adenine dinucleotide phosphate - reduced] + (1.0) nh4_c::[Ammonium] + (1.0) akg_c::[2-Oxoglutarate]
    raw_map["Supp3"]["GLUDy"] = [
        Dict{String, Any}(
            "raw" => [10019.7, 9403.5, 8348.8, 8453.3, 6445.8, 6271.6, 7381.9, 1960.5, 5172.8, 3458.5, 697.4, 501.6], 
            "norm" => [1.279999497, 1.188429919, 1.016800938, 1.034746739, 0.643591514, 0.604065636, 0.839224282, -1.07354627, 0.326177501, -0.254621455, -2.564709642, -3.040158658], 
            "desc" => "MG1655_gdhA_b1761  DB_XREF=PID:g1788059  SEG=NC_000913:+1840395,1841738  LEN=1343  DEF=NADP-specific glutamate dehydrogenase",
        ),
    ]

    # [61] GLUN::[Glutaminase] Glutaminase
    # (-1.0) gln__L_c::[L-Glutamine] + (-1.0) h2o_c::[H2O H2O] ==> (1.0) glu__L_c::[L-Glutamate] + (1.0) nh4_c::[Ammonium]
    raw_map["Supp1"]["GLUN"] = [] 
    raw_map["Supp3"]["GLUN"] = [
        Dict{String, Any}(
            "raw" => [1014.5, 899.9, 1138.6, 1380.4, 1349.4, 1965.9, 1965.2, 4334.2, 3577.5, 2605.3, 7766.7, 7659.8], 
            "norm" => [-1.157057173, -1.32998944, -0.990565034, -0.71273966, -0.745507972, -0.202636101, -0.203149894, 0.937939689, 0.661125728, 0.20362347, 1.779475702, 1.759480685], 
            "desc" => "MG1655_ybaS_b0485  DB_XREF=PID:g1786693  SEG=NC_000913:+510865,511797  LEN=932  DEF=putative glutaminase",
        ),
        Dict{String, Any}(
            "raw" => [838.6, 977.3, 965.1, 758.2, 900.7, 624.0, 703.2, 336.7, 766.2, 445.4, 188.3, 256.8], 
            "norm" => [0.543438373, 0.764257035, 0.74613398, 0.398034, 0.646502204, 0.117001572, 0.289390614, -0.773080736, 0.413176569, -0.369442899, -1.611511457, -1.163899255], 
            "desc" => "MG1655_yneH_b1524  DB_XREF=PID:g1787804  SEG=NC_000913:-1610349,1611275  LEN=926  DEF=putative glutaminase",
        ),
    ]

    # [62] GLUSy::[Glutamate synthase (NADPH)] Glutamate synthase (NADPH)
    # (-1.0) gln__L_c::[L-Glutamine] + (-1.0) h_c::[H+] + (-1.0) nadph_c::[Nicotinamide adenine dinucleotide phosphate - reduced] + (-1.0) akg_c::[2-Oxoglutarate] ==> (2.0) glu__L_c::[L-Glutamate] + (1.0) nadp_c::[Nicotinamide adenine dinucleotide phosphate]
    raw_map["Supp1"]["GLUSy"] = ["gltBD"]
    raw_map["Supp3"]["GLUSy"] = [
        Dict{String, Any}(
            "raw" => [7459.7, 6204.3, 6157.3, 5596.0, 3300.0, 4317.8, 3596.5, 867.4, 1327.2, 643.5, 181.0, 169.2], 
            "norm" => [1.994412095, 1.728562931, 1.717592345, 1.579690446, 0.817760508, 1.205590903, 0.941888087, -1.109936168, -0.496319725, -1.540693463, -3.370643914, -3.467904043], 
            "desc" => "MG1655_gltB_b3212  DB_XREF=PID:g1789605  SEG=NC_000913:+3352267,3356820  LEN=4553  DEF=glutamate synthase, large subunit",
        ),
        Dict{String, Any}(
            "raw" => [4764.0, 3990.5, 3956.6, 3694.2, 1880.7, 3083.8, 2923.3, 1052.7, 992.8, 1653.5, 262.5, 145.8], 
            "norm" => [1.561882802, 1.306278913, 1.29397061, 1.194971363, 0.220979124, 0.934418591, 0.85730728, -0.616196258, -0.700715591, 0.035232434, -2.619901283, -3.468227986], 
            "desc" => "MG1655_gltD_b3213  DB_XREF=PID:g1789606  SEG=NC_000913:+3356833,3358251  LEN=1418  DEF=glutamate synthase, small subunit",
        ),
    ]

    # [63] GLUt2r::[L glutamate transport via proton symport  reversible] L glutamate transport via proton symport  reversible
    # (-1.0) glu__L_e::[L-Glutamate] + (-1.0) h_e::[H+] <==> (1.0) glu__L_c::[L-Glutamate] + (1.0) h_c::[H+]
    raw_map["Supp1"]["GLUt2r"] = []
    raw_map["Supp3"]["GLUt2r"] = []

    # [64] GLYCtpp::[Glycerol transport via channel (periplasm)] Transport, Inner Membrane
    # (-1.0) glyc_e::[Glycerol] <==> (1.0) glyc_c::[Glycerol]
    raw_map["Supp1"]["GLYCtpp"] = []
    raw_map["Supp3"]["GLYCtpp"] = []

    # [65] GLYK::[Glycerol kinase] Alternate Carbon Metabolism
    # (-1.0) atp_c::[ATP C10H12N5O13P3] + (-1.0) glyc_c::[Glycerol] ==> (1.0) h_c::[H+] + (1.0) adp_c::[ADP C10H12N5O10P2] + (1.0) glyc3p_c::[Glycerol 3-phosphate]
    raw_map["Supp1"]["GLYK"] = ["glpK"]
    raw_map["Supp3"]["GLYK"] = [
        Dict{String, Any}(
            "raw" => [28.5, 1.2, 28.0, 36.7, 29.4, 33.2, 10.0, 63.0, 75.5, 44.6, 14.9, 24.4], 
            "norm" => [0.286683571, -4.283172037, 0.261148479, 0.651501715, 0.331537807, 0.506904893, -1.224278348, 1.43107348, 1.692198296, 0.932765362, -0.648966018, 0.0626028], 
            "desc" => "CFT073_glpK_c4878  DB_XREF=26250692  SEG=NC_004431:-4635408,4637021  LEN=1613  DEF=Glycerol kinase",
        ),
        Dict{String, Any}(
            "raw" => [90.5, 72.1, 77.3, 86.0, 1225.7, 1445.7, 265.6, 11978.6, 7582.9, 1490.2, 238.7, 305.4], 
            "norm" => [-2.442828734, -2.770747267, -2.670278112, -2.516409866, 1.316715574, 1.55487787, -0.889563285, 4.605497061, 3.945859361, 1.598615632, -1.043619865, -0.688118369], 
            "desc" => "MG1655_glpK_b3926  DB_XREF=PID:g1790361  SEG=NC_000913:-4113294,4114802  LEN=1508  DEF=glycerol kinase",
        ),
    ]

    # [66] GND::[Phosphogluconate dehydrogenase] Phosphogluconate dehydrogenase
    # (-1.0) nadp_c::[Nicotinamide adenine dinucleotide phosphate] + (-1.0) 6pgc_c::[6-Phospho-D-gluconate] ==> (1.0) nadph_c::[Nicotinamide adenine dinucleotide phosphate - reduced] + (1.0) ru5p__D_c::[D-Ribulose 5-phosphate] + (1.0) co2_c::[CO2 CO2]
    raw_map["Supp1"]["GND"] = ["gnd"]
    raw_map["Supp3"]["GND"] = [
        Dict{String, Any}(
            "raw" => [8693.1, 10265.8, 10241.9, 10632.7, 9287.8, 7664.1, 7102.1, 4456.4, 5131.8, 5372.4, 2991.4, 2597.9], 
            "norm" => [0.442598361, 0.682501774, 0.679139093, 0.733163707, 0.538064526, 0.260844006, 0.150973294, -0.521393647, -0.317807434, -0.251705656, -1.096451544, -1.299926481], 
            "desc" => "MG1655_gnd_b2029  DB_XREF=PID:g1788341  SEG=NC_000913:-2097884,2099290  LEN=1406  DEF=gluconate-6-phosphate dehydrogenase, decarboxylating",
        ),
        Dict{String, Any}(
            "raw" => [11789.5, 13321.7, 13535.5, 14977.9, 12398.6, 10392.6, 10030.4, 6285.5, 6977.3, 6101.4, 3879.8, 3364.7], 
            "norm" => [0.470699729, 0.646975394, 0.669945377, 0.816032558, 0.543374422, 0.288753825, 0.237576335, -0.436703387, -0.286062034, -0.479590585, -1.132748615, -1.338253021], 
            "desc" => "CFT073_gnd_c2556  DB_XREF=26248404  SEG=NC_004431:-2388005,2389411  LEN=1406  DEF=6-phosphogluconate dehydrogenase, decarboxylating",
        ),
        Dict{String, Any}(
            "raw" => [79.1, 77.5, 89.3, 70.2, 96.2, 33.9, 54.9, 87.0, 73.8, 64.0, 49.6, 74.9], 
            "norm" => [0.209498737, 0.180017353, 0.384481218, 0.037292073, 0.491857937, -1.012893684, -0.317372808, 0.346836444, 0.109441859, -0.096107052, -0.463838837, 0.130786761], 
            "desc" => "EDL933_gnd_Z3191  DB_XREF=GI:12516213  SEG=NC_002655:-2842138,2843544  LEN=1406  DEF=gluconate-6-phosphate dehydrogenase, decarboxylating",
        ),
    ]

    # [67] H2Ot::[H2O transport via diffusion] H2O transport via diffusion
    # (-1.0) h2o_e::[H2O H2O] <==> (1.0) h2o_c::[H2O H2O]
    raw_map["Supp1"]["H2Ot"] = [] # Not apply
    raw_map["Supp3"]["H2Ot"] = [] # Not apply

    # [68] ICDHyr::[Isocitrate dehydrogenase (NADP)] Isocitrate dehydrogenase (NADP)
    # (-1.0) icit_c::[Isocitrate] + (-1.0) nadp_c::[Nicotinamide adenine dinucleotide phosphate] <==> (1.0) nadph_c::[Nicotinamide adenine dinucleotide phosphate - reduced] + (1.0) akg_c::[2-Oxoglutarate] + (1.0) co2_c::[CO2 CO2]
    raw_map["Supp1"]["ICDHyr"] = ["icdA"]
    raw_map["Supp3"]["ICDHyr"] = [
        Dict{String, Any}(
            "raw" => [14897.2, 18205.1, 17564.5, 21398.7, 15912.3, 15427.5, 20720.2, 17131.4, 16412.3, 20735.0, 7900.8, 3217.4], 
            "norm" => [0.053059485, 0.342360955, 0.290680799, 0.575541444, 0.148160671, 0.103522585, 0.529056219, 0.254661345, 0.192795721, 0.530086337, -0.861911063, -2.158014496], 
            "desc" => "MG1655_icdA_b1136  DB_XREF=PID:g1787381  SEG=NC_000913:+1194346,1195596  LEN=1250  DEF=isocitrate dehydrogenase, specific for NADP+",
        ),
    ]

    # [69] ICL::[Isocitrate lyase] Isocitrate lyase
    # (-1.0) icit_c::[Isocitrate] ==> (1.0) glx_c::[Glyoxylate] + (1.0) succ_c::[Succinate]
    raw_map["Supp1"]["ICL"] = ["aceA"]
    raw_map["Supp3"]["ICL"] = [
        Dict{String, Any}(
            "raw" => [8530.5, 8403.8, 8347.4, 8465.1, 7557.2, 11577.2, 12365.8, 22640.3, 16599.1, 22203.7, 24274.9, 6526.5], 
            "norm" => [-0.467383157, -0.488971634, -0.498686557, -0.478486351, -0.642161658, -0.026798994, 0.06827021, 0.940807708, 0.493019654, 0.912714739, 1.041379985, -0.853703944], 
            "desc" => "MG1655_aceA_b4015  DB_XREF=PID:g1790445  SEG=NC_000913:+4214688,4215992  LEN=1304  DEF=isocitrate lyase",
        ),
    ]

    # [70] LDH_D::[D-lactate dehydrogenase] D-lactate dehydrogenase
    # (-1.0) lac__D_c::[D-Lactate] + (-1.0) nad_c::[Nicotinamide adenine dinucleotide] <==> (1.0) h_c::[H+] + (1.0) nadh_c::[Nicotinamide adenine dinucleotide - reduced] + (1.0) pyr_c::[Pyruvate]
    raw_map["Supp1"]["LDH_D"] = ["dld"]
    raw_map["Supp3"]["LDH_D"] = [
        Dict{String, Any}(
            "raw" => [5037.4, 4852.6, 3479.9, 3553.6, 4943.6, 4786.6, 4267.2, 3718.0, 3326.8, 3022.8, 1656.8, 947.8], 
            "norm" => [0.601577724, 0.547656374, 0.067944279, 0.098179728, 0.574460446, 0.52789968, 0.36218816, 0.1634252, 0.003033567, -0.135216042, -1.002702111, -1.808447004], 
            "desc" => "MG1655_dld_b2133  DB_XREF=PID:g1788454  SEG=NC_000913:+2220205,2221920  LEN=1715  DEF=D-lactate dehydrogenase, FAD protein, NADH independent",
        ),
        Dict{String, Any}(
            "raw" => [89.0, 61.9, 40.2, 26.8, 63.7, 82.3, 74.5, 70.9, 51.8, 63.9, 68.1, 70.8], 
            "norm" => [0.546425562, 0.022559635, -0.600184273, -1.185146774, 0.063913598, 0.433512656, 0.289860651, 0.218405853, -0.234427676, 0.068436157, 0.160275024, 0.216369586], 
            "desc" => "EDL933_dld_Z3382  DB_XREF=GI:12516441  SEG=NC_002655:+3022365,3024080  LEN=1715  DEF=D-lactate dehydrogenase, FAD protein, NADH independent",
        ),
    ]
    # [71] MALS::[Malate synthase] Malate synthase
    # (-1.0) glx_c::[Glyoxylate] + (-1.0) h2o_c::[H2O H2O] + (-1.0) accoa_c::[Acetyl-CoA] ==> (1.0) h_c::[H+] + (1.0) mal__L_c::[L-Malate] + (1.0) coa_c::[Coenzyme A]
    raw_map["Supp1"]["MALS"] = ["aceB", "glcB"]
    raw_map["Supp3"]["MALS"] = [
        Dict{String, Any}(
            "raw" => [11577.5, 11060.5, 12123.5, 13214.9, 11827.3, 18000.9, 18265.7, 31569.3, 22998.8, 32577.4, 25910.8, 11288.7], 
            "norm" => [-0.547711595, -0.613618746, -0.481229093, -0.356869844, -0.516914587, 0.089033688, 0.110101692, 0.89948692, 0.442523237, 0.944836115, 0.614518208, -0.584155996], 
            "desc" => "MG1655_aceB_b4014  DB_XREF=PID:g1790444  SEG=NC_000913:+4213057,4214658  LEN=1601  DEF=malate synthase A",
        ),
        Dict{String, Any}(
            "raw" => [136.3, 76.3, 89.3, 85.5, 131.7, 248.5, 95.5, 278.3, 182.4, 178.7, 124.3, 80.2], 
            "norm" => [0.072774772, -0.764255828, -0.537278709, -0.600014465, 0.023244556, 0.939235062, -0.440438151, 1.102630119, 0.49309494, 0.463528845, -0.060184493, -0.692336648], 
            "desc" => "CFT073_aceB_c4971  DB_XREF=26250783  SEG=NC_004431:+4746647,4748266  LEN=1619  DEF=Malate synthase A",
        ),
        Dict{String, Any}(
            "raw" => [18.7, 24.6, 14.7, 7.9, 18.4, 14.4, 22.1, 47.3, 26.1, 14.9, 31.9, 3.3], 
            "norm" => [0.13282536, 0.528445405, -0.214396755, -1.110288352, 0.109492856, -0.244144098, 0.37383346, 1.471627273, 0.613836897, -0.194900579, 0.903343514, -2.369674981], 
            "desc" => "CFT073_glcB_c3705  DB_XREF=26249538  SEG=NC_004431:-3536293,3538464  LEN=2171  DEF=Malate synthase G",
        ),
        Dict{String, Any}(
            "raw" => [771.6, 773.3, 709.6, 699.9, 828.2, 6713.6, 2719.3, 5359.5, 2543.7, 1293.5, 6384.8, 1370.2], 
            "norm" => [-1.160160929, -1.156985859, -1.281008062, -1.300865264, -1.058034869, 1.961000606, 0.657149345, 1.636012437, 0.560842554, -0.414805923, 1.888555452, -0.331699487], 
            "desc" => "MG1655_glcB_b2976  DB_XREF=PID:g1789348  SEG=NC_000913:-3119650,3121821  LEN=2171  DEF=malate synthase G",
        ),
    ]

    # [72] MALTabcpp::[Maltose transport via ABC system (periplasm)] Transport, Inner Membrane
    # (-1.0) h2o_c::[H2O H2O] + (-1.0) atp_c::[ATP C10H12N5O13P3] + (-1.0) malt_e::[Maltose C12H22O11] ==> (1.0) h_c::[H+] + (1.0) pi_c::[Phosphate] + (1.0) adp_c::[ADP C10H12N5O10P2] + (1.0) malt_c::[Maltose C12H22O11]
    raw_map["Supp1"]["MALTabcpp"] = []
    raw_map["Supp3"]["MALTabcpp"] = []
    

    # [73] MALt2_2::[Malate transport via proton symport (2 H)] Malate transport via proton symport (2 H)
    # (-2.0) h_e::[H+] + (-1.0) mal__L_e::[L-Malate] ==> (2.0) h_c::[H+] + (1.0) mal__L_c::[L-Malate]
    raw_map["Supp1"]["MALt2_2"] = []
    raw_map["Supp3"]["MALt2_2"] = [
        Dict{String, Any}(
            "raw" => [15.9, 31.8, 15.8, 12.8, 150.5, 134.1, 154.5, 45.3, 1.8, 7.0, 1.2, 3.4], 
            "norm" => [-0.125665351, 0.874334649, -0.134767558, -0.438548306, 3.116999465, 2.950545216, 3.154842817, 1.384818934, -3.268623305, -1.309265289, -3.853585806, -2.351085465], 
            "desc" => "CFT073_malE_c5004  DB_XREF=26250816  SEG=NC_004431:-4779123,4780313  LEN=1190  DEF=Maltose-binding periplasmic protein precursor",
        ),
        Dict{String, Any}(
            "raw" => [1707.8, 2249.7, 1437.2, 1688.0, 20954.5, 19510.7, 24903.7, 11237.2, 323.8, 875.9, 341.6, 466.7], 
            "norm" => [-0.524301664, -0.126708066, -0.773179855, -0.541125791, 3.092747497, 2.989752939, 3.341847502, 2.19377, -2.923265804, -1.487602621, -2.846060815, -2.395873323], 
            "desc" => "MG1655_malE_b4034  DB_XREF=PID:g1790466  SEG=NC_000913:-4242808,4243998  LEN=1190  DEF=periplasmic maltose-binding protein; substrate recognition for transport and chemotaxis",
        ),
        Dict{String, Any}(
            "raw" => [225.9, 193.0, 228.1, 165.4, 5988.0, 6598.9, 4245.4, 551.0, 61.0, 132.2, 75.5, 211.5], 
            "norm" => [-0.84287883, -1.069962253, -0.828896654, -1.292603866, 3.885439216, 4.02559055, 3.389265485, 0.443489218, -2.731681953, -1.615840924, -2.424014551, -0.937905437], 
            "desc" => "MG1655_malG_b4032  DB_XREF=PID:g1790464  SEG=NC_000913:-4240205,4241095  LEN=890  DEF=part of maltose permease, inner membrane",
        ),
        Dict{String, Any}(
            "raw" => [460.2, 593.2, 290.1, 349.6, 9314.4, 8944.5, 6519.7, 648.2, 86.5, 166.7, 80.2, 170.1], 
            "norm" => [-0.367456604, -0.00119899, -1.033167292, -0.764012403, 3.971673346, 3.913211345, 3.457016089, 0.126721433, -2.77894555, -1.832463483, -2.888043446, -1.803334446], 
            "desc" => "MG1655_malF_b4033  DB_XREF=PID:g1790465  SEG=NC_000913:-4241110,4242654  LEN=1544  DEF=part of maltose permease, periplasmic",
        ),
        Dict{String, Any}(
            "raw" => [14.0, 68.6, 29.3, 39.9, 224.6, 131.7, 184.9, 96.8, 13.6, 43.9, 23.1, 50.5], 
            "norm" => [-1.899996379, 0.39278537, -0.834522542, -0.38903446, 2.103862816, 1.333760234, 1.823250113, 0.889583841, -1.941816555, -0.251202267, -1.177530355, -0.049139818], 
            "desc" => "CFT073_malK_c5005  DB_XREF=26250817  SEG=NC_004431:+4780588,4781793  LEN=1205  DEF=Maltosemaltodextrin transport ATP-binding protein malK",
        ),
        Dict{String, Any}(
            "raw" => [473.4, 502.4, 401.9, 418.7, 9472.9, 7983.7, 8217.8, 762.2, 68.1, 347.8, 33.4, 155.1], 
            "norm" => [-0.349883858, -0.2641071, -0.586106987, -0.527026646, 3.972790685, 3.726042041, 3.76773675, 0.337226044, -3.147216861, -0.794685632, -4.175023557, -1.959744878], 
            "desc" => "MG1655_malK_b4035  DB_XREF=PID:g1790467  SEG=NC_000913:+4244363,4245478  LEN=1115  DEF=ATP-binding component of transport system for maltose",
        ),
    ]

    # [74] MDH::[Malate dehydrogenase] Malate dehydrogenase
    # (-1.0) mal__L_c::[L-Malate] + (-1.0) nad_c::[Nicotinamide adenine dinucleotide] <==> (1.0) h_c::[H+] + (1.0) nadh_c::[Nicotinamide adenine dinucleotide - reduced] + (1.0) oaa_c::[Oxaloacetate]
    raw_map["Supp1"]["MDH"] = ["mdh"]
    raw_map["Supp3"]["MDH"] = [
        Dict{String, Any}(
            "raw" => [12702.1, 14957.3, 15624.4, 19966.1, 14667.8, 16270.0, 20483.4, 21576.8, 18494.2, 25217.5, 14618.1, 8209.4], 
            "norm" => [-0.35915594, -0.123373201, -0.060422184, 0.293329584, -0.151570473, -0.002008722, 0.330232232, 0.405257945, 0.182849922, 0.630202285, -0.156467165, -0.988874284], 
            "desc" => "MG1655_mdh_b3236  DB_XREF=PID:g1789632  SEG=NC_000913:-3380965,3381903  LEN=938  DEF=malate dehydrogenase",
        ),
    ]

    # [75] ME1::[Malic enzyme (NAD)] Malic enzyme (NAD)
    # (-1.0) mal__L_c::[L-Malate] + (-1.0) nad_c::[Nicotinamide adenine dinucleotide] ==> (1.0) nadh_c::[Nicotinamide adenine dinucleotide - reduced] + (1.0) pyr_c::[Pyruvate] + (1.0) co2_c::[CO2 CO2]
    raw_map["Supp1"]["ME1"] = ["maeA"]
    raw_map["Supp3"]["ME1"] = [
        Dict{String, Any}(
            "raw" => [1996.6, 2151.3, 2068.3, 1913.2, 2294.7, 1450.7, 1435.2, 759.0, 1327.5, 376.6, 307.3, 253.1], 
            "norm" => [0.892463548, 1.00012694, 0.943363675, 0.830905913, 1.09322377, 0.431667422, 0.416170012, -0.502909993, 0.303630077, -1.513976879, -1.807362112, -2.087302371], 
            "desc" => "MG1655_sfcA_b1479  DB_XREF=PID:g1787754  SEG=NC_000913:-1551996,1553720  LEN=1724  DEF=NAD-linked malate dehydrogenase (malic enzyme)",
        ),
    ]

    # [76] ME2::[Malic enzyme (NADP)] Malic enzyme (NADP)
    # (-1.0) mal__L_c::[L-Malate] + (-1.0) nadp_c::[Nicotinamide adenine dinucleotide phosphate] ==> (1.0) nadph_c::[Nicotinamide adenine dinucleotide phosphate - reduced] + (1.0) pyr_c::[Pyruvate] + (1.0) co2_c::[CO2 CO2]
    raw_map["Supp1"]["ME2"] = ["maeB"]
    raw_map["Supp3"]["ME2"] = [
        Dict{String, Any}(
            "raw" => [41.7, 37.2, 69.6, 71.0, 40.9, 43.2, 95.5, 121.0, 71.7, 143.4, 143.2, 130.4], 
            "norm" => [-0.840502796, -1.005247558, -0.101462874, -0.072731155, -0.868449337, -0.789518867, 0.354950553, 0.696384963, -0.058577061, 0.941422939, 0.939409408, 0.804321785], 
            "desc" => "CFT073_c2988  DB_XREF=26248830  SEG=NC_004431:-2846856,2849135  LEN=2279  DEF=NADP-dependent malic enzyme",
        ),
    ]

    # [77] NADH16::[NADH dehydrogenase (ubiquinone-8 & 3 protons)] NADH dehydrogenase (ubiquinone-8 & 3 protons)
    # (-4.0) h_c::[H+] + (-1.0) nadh_c::[Nicotinamide adenine dinucleotide - reduced] + (-1.0) q8_c::[Ubiquinone-8] ==> (3.0) h_e::[H+] + (1.0) nad_c::[Nicotinamide adenine dinucleotide] + (1.0) q8h2_c::[Ubiquinol-8]
    raw_map["Supp1"]["NADH16"] = ["nuoABEFGHIJKLMN"]
    raw_map["Supp3"]["NADH16"] = [
        Dict{String, Any}(
            "raw" => [5223.6, 5515.1, 4997.8, 5658.3, 7376.8, 6653.7, 5647.1, 5592.4, 4669.9, 6367.4, 4356.5, 3365.3], 
            "norm" => [-0.034640883, 0.043701737, -0.098392139, 0.080683361, 0.463309813, 0.314471511, 0.077824871, 0.063782246, -0.196293652, 0.251019089, -0.296515765, -0.668950189], 
            "desc" => "MG1655_nuoA_b2288  DB_XREF=PID:g1788625  SEG=NC_000913:-2402649,2403092  LEN=443  DEF=NADH dehydrogenase I chain A",
        ),
        Dict{String, Any}(
            "raw" => [4779.2, 5061.4, 5028.2, 5635.0, 6350.6, 6207.0, 4675.0, 4212.7, 3771.2, 4725.9, 2684.2, 2191.1], 
            "norm" => [0.112350464, 0.195117816, 0.185623356, 0.349996932, 0.522474224, 0.489477467, 0.080547686, -0.069673498, -0.229395015, 0.096170424, -0.719926508, -1.012763349], 
            "desc" => "MG1655_nuoB_b2287  DB_XREF=PID:g1788624  SEG=NC_000913:-2401971,2402633  LEN=662  DEF=NADH dehydrogenase I chain B",
        ),
        Dict{String, Any}(
            "raw" => [5722.4, 5464.6, 5957.1, 6645.4, 7139.5, 6238.7, 5620.5, 5735.9, 4094.0, 4614.9, 3370.5, 1932.4], 
            "norm" => [0.207899, 0.14139455, 0.265888831, 0.423634694, 0.527101695, 0.332524089, 0.181977132, 0.211298525, -0.275210244, -0.10242196, -0.555758722, -1.358327589], 
            "desc" => "MG1655_nuoE_b2285  DB_XREF=PID:g1788621  SEG=NC_000913:-2399572,2400072  LEN=500  DEF=NADH dehydrogenase I chain E",
        ),
        Dict{String, Any}(
            "raw" => [65.0, 56.2, 106.2, 68.8, 84.1, 90.6, 164.6, 65.5, 68.5, 119.1, 75.0, 71.7], 
            "norm" => [-0.339993921, -0.549863509, 0.368278222, -0.258025075, 0.031672161, 0.139077411, 1.000458791, -0.328938733, -0.264329651, 0.533667869, -0.133543044, -0.198460521], 
            "desc" => "EDL933_nuoF_Z3543  DB_XREF=GI:12516634  SEG=NC_002655:-3196344,3197681  LEN=1337  DEF=NADH dehydrogenase I chain F",
        ),
        Dict{String, Any}(
            "raw" => [3898.7, 3080.6, 3473.8, 3160.4, 4162.1, 3836.9, 3384.0, 2331.7, 2209.5, 1645.9, 1732.6, 879.3], 
            "norm" => [0.588813685, 0.249031907, 0.422335234, 0.285927706, 0.683132167, 0.565761704, 0.384550108, -0.15279728, -0.23045953, -0.655302776, -0.581240838, -1.559752087], 
            "desc" => "MG1655_nuoF_b2284  DB_XREF=PID:g1788620  SEG=NC_000913:-2398238,2399575  LEN=1337  DEF=NADH dehydrogenase I chain F",
        ),
        Dict{String, Any}(
            "raw" => [26.7, 49.6, 44.4, 33.3, 32.3, 44.5, 57.4, 68.3, 37.4, 50.8, 56.0, 82.4], 
            "norm" => [-0.793783213, 0.099717165, -0.060063279, -0.475100778, -0.51908879, -0.056817619, 0.310427782, 0.561262623, -0.307584685, 0.134205542, 0.274803872, 0.832021382], 
            "desc" => "CFT073_nuoG_c2824  DB_XREF=26248670  SEG=NC_004431:-2679944,2682676  LEN=2732  DEF=NADH dehydrogenase I chain G",
        ),
        Dict{String, Any}(
            "raw" => [1338.1, 1121.0, 1194.5, 941.2, 1408.3, 1329.0, 875.7, 540.7, 595.5, 383.7, 481.4, 335.3], 
            "norm" => [0.768626479, 0.513226821, 0.604847396, 0.261013768, 0.842395236, 0.758781647, 0.156949159, -0.538659196, -0.399386044, -1.033508787, -0.706251411, -1.228035069], 
            "desc" => "MG1655_nuoG_b2283  DB_XREF=PID:g1788619  SEG=NC_000913:-2395459,2398191  LEN=2732  DEF=NADH dehydrogenase I chain G",
        ),
        Dict{String, Any}(
            "raw" => [55.0, 37.6, 47.9, 28.4, 30.7, 33.6, 24.2, 33.1, 20.7, 14.6, 22.3, 25.6], 
            "norm" => [0.908314479, 0.359615522, 0.708908516, -0.04522621, 0.067121516, 0.197344093, -0.276110092, 0.175714077, -0.501486372, -1.005148771, -0.39407343, -0.194973329], 
            "desc" => "EDL933_nuoG_Z3542  DB_XREF=GI:12516633  SEG=NC_002655:-3193565,3196297  LEN=2732  DEF=NADH dehydrogenase I chain G",
        ),
        Dict{String, Any}(
            "raw" => [2314.1, 2113.1, 2271.8, 2077.7, 3110.3, 2490.3, 1741.2, 1203.8, 859.1, 622.7, 851.7, 426.3], 
            "norm" => [0.694988563, 0.563898396, 0.668373185, 0.539524712, 1.121591094, 0.800856904, 0.284619279, -0.247866925, -0.734564669, -1.198853462, -0.747045391, -1.745521686], 
            "desc" => "MG1655_nuoH_b2282  DB_XREF=PID:g1788618  SEG=NC_000913:-2394485,2395462  LEN=977  DEF=NADH dehydrogenase I chain H",
        ),
        Dict{String, Any}(
            "raw" => [3796.0, 3448.5, 3647.1, 3231.3, 4229.9, 3575.2, 3087.0, 1896.3, 1629.5, 1050.9, 1330.1, 730.3], 
            "norm" => [0.718521307, 0.580010282, 0.660791072, 0.486156014, 0.874664871, 0.632065265, 0.420246798, -0.282771465, -0.501529333, -1.134333291, -0.79442397, -1.65939755], 
            "desc" => "MG1655_nuoI_b2281  DB_XREF=PID:g1788617  SEG=NC_000913:-2393928,2394470  LEN=542  DEF=NADH dehydrogenase I chain I",
        ),
        Dict{String, Any}(
            "raw" => [5071.3, 4856.7, 4601.8, 5737.9, 6247.7, 5424.0, 4216.7, 3077.6, 2605.5, 2128.6, 2222.2, 1228.5], 
            "norm" => [0.505188353, 0.442809104, 0.365031015, 0.683355556, 0.806157911, 0.602189909, 0.238947115, -0.215361534, -0.455607014, -0.747262401, -0.685178603, -1.540269411], 
            "desc" => "MG1655_nuoJ_b2280  DB_XREF=PID:g1788616  SEG=NC_000913:-2393362,2393916  LEN=554  DEF=NADH dehydrogenase I chain J",
        ),
        Dict{String, Any}(
            "raw" => [3536.2, 3340.0, 3366.9, 3204.3, 4343.0, 3091.2, 2915.9, 1738.9, 1624.0, 1214.4, 1173.8, 600.6], 
            "norm" => [0.69237391, 0.610022139, 0.621594909, 0.550183262, 0.988865989, 0.498341036, 0.41411528, -0.331650994, -0.430274331, -0.849582267, -0.89863935, -1.865349583], 
            "desc" => "MG1655_nuoK_b2279  DB_XREF=PID:g1788615  SEG=NC_000913:-2393063,2393365  LEN=302  DEF=NADH dehydrogenase I chain K",
        ),
        Dict{String, Any}(
            "raw" => [4689.1, 4245.1, 4254.8, 3702.7, 4895.5, 3846.1, 3620.2, 1816.8, 2114.8, 1402.7, 1280.7, 780.6], 
            "norm" => [0.828193069, 0.684680562, 0.687973341, 0.487459687, 0.890338239, 0.542278296, 0.454951425, -0.539718366, -0.320596745, -0.912911489, -1.044185409, -1.758462609], 
            "desc" => "MG1655_nuoL_b2278  DB_XREF=PID:g1788614  SEG=NC_000913:-2391225,2393066  LEN=1841  DEF=NADH dehydrogenase I chain L",
        ),
        Dict{String, Any}(
            "raw" => [635.8, 449.1, 552.9, 414.5, 829.0, 547.7, 425.5, 268.3, 324.3, 178.0, 180.5, 103.2], 
            "norm" => [0.852038629, 0.350502334, 0.650484184, 0.234837714, 1.234837714, 0.636851493, 0.272624744, -0.392687333, -0.119205364, -0.984657147, -0.964535551, -1.771091417], 
            "desc" => "MG1655_nuoM_b2277  DB_XREF=PID:g1788613  SEG=NC_000913:-2389532,2391061  LEN=1529  DEF=NADH dehydrogenase I chain M",
        ),
        Dict{String, Any}(
            "raw" => [2638.4, 2646.0, 2688.6, 2401.1, 2858.9, 2007.9, 3177.8, 1413.6, 1690.2, 1266.0, 1459.1, 1051.6], 
            "norm" => [0.408187219, 0.412336976, 0.435379047, 0.272219404, 0.523984072, 0.014211335, 0.676552244, -0.49210214, -0.234282116, -0.65119868, -0.446397323, -0.918890038], 
            "desc" => "MG1655_nuoN_b2276  DB_XREF=PID:g1788612  SEG=NC_000913:-2388068,2389345  LEN=1277  DEF=NADH dehydrogenase I chain N",
        ),
    ]

    # [78] NADTRHD::[NAD transhydrogenase] NAD transhydrogenase
    # (-1.0) nad_c::[Nicotinamide adenine dinucleotide] + (-1.0) nadph_c::[Nicotinamide adenine dinucleotide phosphate - reduced] ==> (1.0) nadh_c::[Nicotinamide adenine dinucleotide - reduced] + (1.0) nadp_c::[Nicotinamide adenine dinucleotide phosphate]
    raw_map["Supp1"]["NADTRHD"] = ["pntAB"]
    raw_map["Supp3"]["NADTRHD"] = [
        Dict{String, Any}(
            "raw" => [13071.5, 14675.0, 14760.1, 18214.6, 10825.3, 10443.5, 11896.1, 5841.4, 9354.4, 8085.3, 3590.7, 2353.8], 
            "norm" => [0.548844622, 0.715780421, 0.724122413, 1.027515231, 0.276826924, 0.22502521, 0.412908597, -0.613193998, 0.066136943, -0.144206872, -1.315243055, -1.924516436], 
            "desc" => "MG1655_pntA_b1603  DB_XREF=PID:g1787887  SEG=NC_000913:-1674395,1675927  LEN=1532  DEF=pyridine nucleotide transhydrogenase, alpha subunit",
        ),
        Dict{String, Any}(
            "raw" => [9321.5, 9529.0, 10441.3, 10384.0, 6494.6, 6282.7, 6642.8, 1982.1, 4586.7, 3669.2, 1497.2, 993.0], 
            "norm" => [0.974903126, 1.006665818, 1.138570438, 1.13063138, 0.45358167, 0.405725689, 0.486132475, -1.258629253, -0.048202455, -0.370193459, -1.66339205, -2.255793381], 
            "desc" => "MG1655_pntB_b1602  DB_XREF=PID:g1787886  SEG=NC_000913:-1672996,1674384  LEN=1388  DEF=pyridine nucleotide transhydrogenase, beta subunit",
        ),
    ]

    # [79] NH4t::[Ammonia reversible transport] Ammonia reversible transport
    # (-1.0) nh4_e::[Ammonium] <==> (1.0) nh4_c::[Ammonium]
    raw_map["Supp1"]["NH4t"] = [] # Not found
    raw_map["Supp3"]["NH4t"] = [] # Not found


    # [80] O2t::[O2 transport diffusion ] O2 transport  diffusion 
    # (-1.0) o2_e::[O2 O2] <==> (1.0) o2_c::[O2 O2]
    raw_map["Supp1"]["O2t"] = [] # Not apply
    raw_map["Supp3"]["O2t"] = [] # Not apply

    # [81] PDH::[Pyruvate dehydrogenase] Pyruvate dehydrogenase
    # (-1.0) nad_c::[Nicotinamide adenine dinucleotide] + (-1.0) pyr_c::[Pyruvate] + (-1.0) coa_c::[Coenzyme A] ==> (1.0) nadh_c::[Nicotinamide adenine dinucleotide - reduced] + (1.0) accoa_c::[Acetyl-CoA] + (1.0) co2_c::[CO2 CO2]
    raw_map["Supp1"]["PDH"] = ["lpdA", "aceEF"]
    raw_map["Supp3"]["PDH"] = [
        Dict{String, Any}(
            "raw" => [12674.3, 13320.1, 12464.3, 13417.8, 17059.9, 5347.1, 3858.4, 1594.7, 3032.5, 2692.4, 1333.9, 809.0], 
            "norm" => [1.409082698, 1.480781541, 1.384978491, 1.491324773, 1.837785818, 0.164005191, -0.306748752, -1.581466422, -0.654243822, -0.825858705, -1.839100952, -2.560539859], 
            "desc" => "MG1655_aceE_b0114  DB_XREF=PID:g1786304  SEG=NC_000913:+123017,125680  LEN=2663  DEF=pyruvate dehydrogenase (decarboxylase component)",
        ),
        Dict{String, Any}(
            "raw" => [8927.6, 9869.7, 8826.1, 9282.6, 12365.6, 5235.5, 2893.8, 1622.6, 2427.6, 4986.3, 1386.8, 810.9], 
            "norm" => [1.088914055, 1.233647899, 1.072417759, 1.145170618, 1.558902004, 0.31896899, -0.536393118, -1.371050941, -0.789827608, 0.248611351, -1.597598592, -2.371762416], 
            "desc" => "EDL933_aceF_Z0125  DB_XREF=GI:12512824  SEG=NC_002655:+130185,132077  LEN=1892  DEF=pyruvate dehydrogenase (dihydrolipoyltransacetylase component)",
        ),
        Dict{String, Any}(
            "raw" => [13646.1, 16078.3, 14501.5, 17667.9, 16382.5, 13625.9, 12624.3, 6149.4, 7986.4, 8169.9, 3066.4, 1511.2], 
            "norm" => [0.595442346, 0.832068527, 0.68315579, 0.968084225, 0.859109184, 0.593305177, 0.483157047, -0.554528789, -0.177429111, -0.144656022, -1.558428539, -2.579279835], 
            "desc" => "MG1655_lpdA_b0116  DB_XREF=PID:g1786307  SEG=NC_000913:+127912,129336  LEN=1424  DEF=lipoamide dehydrogenase (NADH); component of 2-oxodehydrogenase and pyruvate complexes; L-protein of glycine cleavage complex",
        ),
    ]

    # [82] PFK::[Phosphofructokinase] Phosphofructokinase
    # (-1.0) atp_c::[ATP C10H12N5O13P3] + (-1.0) f6p_c::[D-Fructose 6-phosphate] ==> (1.0) h_c::[H+] + (1.0) adp_c::[ADP C10H12N5O10P2] + (1.0) fdp_c::[D-Fructose 1,6-bisphosphate]
    raw_map["Supp1"]["PFK"] = ["pfkA", "pfkB"]
    raw_map["Supp3"]["PFK"] = [
        Dict{String, Any}(
            "raw" => [4490.7, 5251.7, 4745.8, 5038.4, 5800.3, 7433.8, 7715.4, 9062.1, 4375.4, 5924.4, 10269.3, 10995.3], 
            "norm" => [-0.517419699, -0.29157554, -0.437708742, -0.351394383, -0.148232525, 0.209739829, 0.263380909, 0.495485366, -0.554945131, -0.117690995, 0.675905894, 0.774455017], 
            "desc" => "CFT073_pfkA_c4867  DB_XREF=26250681  SEG=NC_004431:+4627114,4628172  LEN=1058  DEF=6-phosphofructokinase isozyme I",
        ),
        Dict{String, Any}(
            "raw" => [4094.0, 4304.5, 4584.7, 5455.4, 4706.9, 6353.5, 6556.7, 7904.9, 3507.4, 3715.4, 7100.0, 7177.3], 
            "norm" => [-0.362913456, -0.290578891, -0.199597225, 0.051260423, -0.161647356, 0.271127, 0.316545327, 0.586322653, -0.586022587, -0.502907019, 0.431394466, 0.447016666], 
            "desc" => "MG1655_pfkA_b3916  DB_XREF=PID:g1790350  SEG=NC_000913:+4105132,4106094  LEN=962  DEF=6-phosphofructokinase I",
        ),
        Dict{String, Any}(
            "raw" => [3274.0, 3198.3, 3385.7, 3459.8, 3095.2, 4105.9, 3270.0, 3945.9, 2922.8, 1962.5, 3204.9, 3170.3], 
            "norm" => [0.030985904, -0.002763148, 0.079385727, 0.110620225, -0.050035784, 0.357630073, 0.029222218, 0.300285977, -0.132717305, -0.707375763, 0.000210925, -0.015449051], 
            "desc" => "MG1655_pfkB_b1723  DB_XREF=PID:g1788017  SEG=NC_000913:+1804394,1805323  LEN=929  DEF=6-phosphofructokinase II; suppressor of pfkA",
        ),
    ]

    # [83] PFL::[Pyruvate formate lyase] Pyruvate formate lyase
    # (-1.0) pyr_c::[Pyruvate] + (-1.0) coa_c::[Coenzyme A] ==> (1.0) accoa_c::[Acetyl-CoA] + (1.0) for_c::[Formate]
    raw_map["Supp1"]["PFL"] = ["pflAB", "pflCD"]
    raw_map["Supp3"]["PFL"] = [
        Dict{String, Any}(
            "raw" => [1436.5, 1779.8, 1630.5, 1780.3, 1752.5, 1832.9, 1848.8, 2096.8, 1127.4, 1625.1, 2426.9, 3269.1], 
            "norm" => [-0.343405759, -0.03424862, -0.160649311, -0.03384338, -0.056549308, 0.008164325, 0.020625413, 0.202225507, -0.692964279, -0.165435255, 0.413150912, 0.842929757], 
            "desc" => "MG1655_pflA_b0902  DB_XREF=PID:g1787130  SEG=NC_000913:-949563,950303  LEN=740  DEF=pyruvate formate lyase activating enzyme 1",
        ),
        Dict{String, Any}(
            "raw" => [10198.4, 11748.9, 10908.1, 14596.2, 11567.3, 9578.9, 11506.1, 15905.6, 10766.9, 14414.8, 18121.7, 20581.6], 
            "norm" => [-0.345407681, -0.141224821, -0.248350679, 0.171842314, -0.163698356, -0.435818613, -0.171351597, 0.295784284, -0.267147581, 0.153800309, 0.483967791, 0.66760463], 
            "desc" => "EDL933_pflB_Z1248  DB_XREF=GI:12514066  SEG=NC_002655:-1170138,1172420  LEN=2282  DEF=formate acetyltransferase 1",
        ),
        Dict{String, Any}(
            "raw" => [33.4, 35.4, 59.4, 29.0, 31.8, 36.4, 49.4, 36.3, 30.9, 82.2, 43.8, 72.9], 
            "norm" => [-0.346081399, -0.262180141, 0.48453343, -0.549876601, -0.416902736, -0.221991051, 0.21858154, -0.225959953, -0.458322663, 0.953208893, 0.045001368, 0.779989313], 
            "desc" => "MG1655_pflC_b3952  DB_XREF=PID:g1790389  SEG=NC_000913:+4143837,4144715  LEN=878  DEF=probable pyruvate formate lyase activating enzyme 2",
        ),
        Dict{String, Any}(
            "raw" => [77.0, 64.3, 99.9, 64.5, 70.3, 95.3, 93.9, 88.5, 65.2, 156.3, 128.3, 182.4], 
            "norm" => [-0.27379701, -0.533836718, 0.101829223, -0.529356295, -0.405130766, 0.033820759, 0.012469702, -0.072978, -0.513783491, 0.747590418, 0.46279381, 0.970378369], 
            "desc" => "MG1655_pflD_b3951  DB_XREF=PID:g1790388  SEG=NC_000913:+4141574,4143871  LEN=2297  DEF=formate acetyltransferase 2",
        ),
    ]

    # [84] PGI::[Glucose-6-phosphate isomerase] Glucose-6-phosphate isomerase
    # (-1.0) g6p_c::[D-Glucose 6-phosphate] <==> (1.0) f6p_c::[D-Fructose 6-phosphate]
    raw_map["Supp1"]["PGI"] = ["pgi"]
    raw_map["Supp3"]["PGI"] = [
        Dict{String, Any}(
            "raw" => [5671.6, 5761.6, 4867.0, 5944.6, 4533.4, 4980.5, 4580.6, 4512.1, 3670.1, 1611.6, 2423.5, 2691.1], 
            "norm" => [0.502034234, 0.524747951, 0.281311221, 0.569858184, 0.178871908, 0.31456903, 0.193815032, 0.172077488, -0.125902181, -1.313227843, -0.724629471, -0.573525552], 
            "desc" => "MG1655_pgi_b4025  DB_XREF=PID:g1790457  SEG=NC_000913:+4231337,4232986  LEN=1649  DEF=glucosephosphate isomerase",
        ),
    ]

    # [85] PGK::[Phosphoglycerate kinase] Phosphoglycerate kinase
    # (-1.0) 3pg_c::[3-Phospho-D-glycerate] + (-1.0) atp_c::[ATP C10H12N5O13P3] <==> (1.0) 13dpg_c::[3-Phospho-D-glyceroyl phosphate] + (1.0) adp_c::[ADP C10H12N5O10P2]
    raw_map["Supp1"]["PGK"] = ["pgk"]
    raw_map["Supp3"]["PGK"] = [
        Dict{String, Any}(
            "raw" => [8712.2, 8728.7, 9095.9, 8901.4, 6974.4, 8807.3, 8933.8, 8108.2, 7712.6, 4239.7, 7357.3, 7436.9], 
            "norm" => [0.163659211, 0.166388941, 0.225838531, 0.194654397, -0.157308753, 0.179321947, 0.199896095, 0.060003813, -0.012160572, -0.875415678, -0.080201443, -0.064676489], 
            "desc" => "MG1655_pgk_b2926  DB_XREF=PID:g1789294  SEG=NC_000913:-3069479,3070642  LEN=1163  DEF=phosphoglycerate kinase",
        ),
    ]

    # [86] PGL::[6-phosphogluconolactonase] 6-phosphogluconolactonase
    # (-1.0) h2o_c::[H2O H2O] + (-1.0) 6pgl_c::[6-phospho-D-glucono-1,5-lactone] ==> (1.0) h_c::[H+] + (1.0) 6pgc_c::[6-Phospho-D-gluconate]
    raw_map["Supp1"]["PGL"] = ["pgl"] 
    raw_map["Supp3"]["PGL"] = [] # Not found 

    # [87] PGM::[Phosphoglycerate mutase] Phosphoglycerate mutase
    # (-1.0) 2pg_c::[D-Glycerate 2-phosphate] <==> (1.0) 3pg_c::[3-Phospho-D-glycerate]
    raw_map["Supp1"]["PGM"] = ["gpmA", "gpmB"]
    raw_map["Supp3"]["PGM"] = [
        Dict{String, Any}(
            "raw" => [12614.6, 15076.4, 13629.1, 18113.0, 11945.1, 12842.5, 20184.4, 24396.9, 21393.0, 30269.1, 29214.6, 34034.9], 
            "norm" => [-0.594661773, -0.337464256, -0.483065936, -0.072730718, -0.673337301, -0.568830159, 0.083484469, 0.35694161, 0.167382577, 0.668089547, 0.616933303, 0.837258638], 
            "desc" => "MG1655_gpmA_b0755  DB_XREF=PID:g1786970  SEG=NC_000913:-786066,786818  LEN=752  DEF=phosphoglyceromutase 1",
        ),
        Dict{String, Any}(
            "raw" => [736.2, 681.6, 639.9, 662.8, 505.5, 411.7, 557.9, 298.2, 600.7, 421.7, 328.3, 413.5], 
            "norm" => [0.553697894, 0.44252548, 0.351446611, 0.402173747, 0.011311236, -0.284806407, 0.153606696, -0.750119598, 0.260244808, -0.250182834, -0.611385106, -0.278512526], 
            "desc" => "MG1655_gpmB_b4395  DB_XREF=PID:g1790856  SEG=NC_000913:+4631366,4632013  LEN=647  DEF=phosphoglyceromutase 2",
        ),
    ]

    # [88] PGMT::[Phosphoglucomutase] Alternate Carbon Metabolism
    # (-1.0) g1p_c::[D-Glucose 1-phosphate] <==> (1.0) g6p_c::[D-Glucose 6-phosphate]
    raw_map["Supp1"]["PGMT"] = ["pgm"]
    raw_map["Supp3"]["PGMT"] = [
        Dict{String, Any}(
            "raw" => [138.9, 144.2, 196.3, 199.1, 129.0, 211.7, 233.8, 189.8, 207.0, 286.7, 240.5, 212.1], 
            "norm" => [-0.483267784, -0.429243219, 0.015745789, 0.036178838, -0.589943318, 0.124706886, 0.267960546, -0.032834391, 0.092316384, 0.562227521, 0.30872251, 0.127430237], 
            "desc" => "CFT073_pgm_c0775  DB_XREF=26246665  SEG=NC_004431:+755903,757594  LEN=1691  DEF=Phosphoglucomutase",
        ),
        Dict{String, Any}(
            "raw" => [2829.6, 2759.3, 2976.0, 2637.2, 2860.0, 3097.7, 2771.4, 2676.6, 2793.8, 2751.8, 1850.8, 1173.0], 
            "norm" => [0.164225839, 0.127930034, 0.237002241, 0.062634701, 0.179642862, 0.294825146, 0.134242667, 0.084029269, 0.145856461, 0.124003334, -0.448223281, -1.106169272], 
            "desc" => "MG1655_pgm_b0688  DB_XREF=PID:g1786904  SEG=NC_000913:+712781,714421  LEN=1640  DEF=phosphoglucomutase",
        ),
    ]

    # [89] PIt2r::[Phosphate reversible transport via symport] Phosphate reversible transport via symport
    # (-1.0) h_e::[H+] + (-1.0) pi_e::[Phosphate] <==> (1.0) h_c::[H+] + (1.0) pi_c::[Phosphate]
    raw_map["Supp1"]["PIt2r"] = [] # Not found
    raw_map["Supp3"]["PIt2r"] = [] # Not found

    # [90] PPC::[Phosphoenolpyruvate carboxylase] Phosphoenolpyruvate carboxylase
    # (-1.0) h2o_c::[H2O H2O] + (-1.0) pep_c::[Phosphoenolpyruvate] + (-1.0) co2_c::[CO2 CO2] ==> (1.0) h_c::[H+] + (1.0) oaa_c::[Oxaloacetate] + (1.0) pi_c::[Phosphate]
    raw_map["Supp1"]["PPC"] = ["ppc"]
    raw_map["Supp3"]["PPC"] = [
        Dict{String, Any}(
            "raw" => [8449.7, 8390.6, 8023.6, 7906.6, 5499.4, 5680.2, 5809.3, 2100.9, 2946.6, 1713.3, 610.7, 711.5], 
            "norm" => [1.258661001, 1.24853486, 1.184010565, 1.16281832, 0.639035105, 0.685702608, 0.718125215, -0.749231626, -0.261187891, -1.043461329, -2.531703369, -2.311303458], 
            "desc" => "MG1655_ppc_b3956  DB_XREF=PID:g1790393  SEG=NC_000913:-4148026,4150677  LEN=2651  DEF=phosphoenolpyruvate carboxylase",
        ),
    ]

    # [91] PPCK::[Phosphoenolpyruvate carboxykinase] Phosphoenolpyruvate carboxykinase
    # (-1.0) oaa_c::[Oxaloacetate] + (-1.0) atp_c::[ATP C10H12N5O13P3] ==> (1.0) pep_c::[Phosphoenolpyruvate] + (1.0) adp_c::[ADP C10H12N5O10P2] + (1.0) co2_c::[CO2 CO2]
    raw_map["Supp1"]["PPCK"] = ["pckA"]
    raw_map["Supp3"]["PPCK"] = [
        Dict{String, Any}(
            "raw" => [1798.6, 1722.1, 1912.7, 1980.8, 7143.3, 7234.4, 3955.3, 8423.5, 6557.4, 13008.8, 5702.4, 1979.9], 
            "norm" => [-1.170954361, -1.233659815, -1.082218125, -1.031745516, 0.818761977, 0.837044634, -0.034041611, 1.056591069, 0.695295166, 1.683587246, 0.493740507, -1.03240117], 
            "desc" => "MG1655_pckA_b3403  DB_XREF=PID:g1789807  SEG=NC_000913:+3530456,3532078  LEN=1622  DEF=phosphoenolpyruvate carboxykinase",
        ),
    ]

    # [92] PPS::[Phosphoenolpyruvate synthase] Phosphoenolpyruvate synthase
    # (-1.0) h2o_c::[H2O H2O] + (-1.0) pyr_c::[Pyruvate] + (-1.0) atp_c::[ATP C10H12N5O13P3] ==> (2.0) h_c::[H+] + (1.0) pep_c::[Phosphoenolpyruvate] + (1.0) pi_c::[Phosphate] + (1.0) amp_c::[AMP C10H12N5O7P]
    raw_map["Supp1"]["PPS"] = ["ppsA"]
    raw_map["Supp3"]["PPS"] = [
        Dict{String, Any}(
            "raw" => [5717.0, 6150.9, 7300.4, 5627.5, 6955.4, 2641.8, 1773.5, 2193.7, 5414.1, 4127.0, 2729.7, 1617.2], 
            "norm" => [0.575241518, 0.680780748, 0.927958741, 0.552477379, 0.858116713, -0.538495523, -1.113417442, -0.80665053, 0.496704762, 0.105076666, -0.491274368, -1.246518664], 
            "desc" => "MG1655_ppsA_b1702  DB_XREF=PID:g1787994  SEG=NC_000913:-1782758,1785136  LEN=2378  DEF=phosphoenolpyruvate synthase",
        ),
    ]

    # [93] PTAr::[Phosphotransacetylase] Phosphotransacetylase
    # (-1.0) pi_c::[Phosphate] + (-1.0) accoa_c::[Acetyl-CoA] <==> (1.0) actp_c::[Acetyl phosphate] + (1.0) coa_c::[Coenzyme A]
    raw_map["Supp1"]["PTAr"] = ["pta"]
    raw_map["Supp3"]["PTAr"] = [
        Dict{String, Any}(
            "raw" => [3295.5, 2969.3, 3065.9, 3072.3, 6874.0, 2099.3, 1486.6, 1774.0, 2425.2, 1318.5, 2154.9, 1364.1], 
            "norm" => [0.46975551, 0.319381001, 0.365568781, 0.368577237, 1.530407991, -0.180833511, -0.678725347, -0.423735851, 0.027361867, -0.851844289, -0.143120939, -0.802792451], 
            "desc" => "MG1655_pta_b2297  DB_XREF=PID:g1788635  SEG=NC_000913:+2412767,2414911  LEN=2144  DEF=phosphotransacetylase",
        ),
    ]

    # [94] PYK::[Pyruvate kinase] Pyruvate kinase
    # (-1.0) h_c::[H+] + (-1.0) pep_c::[Phosphoenolpyruvate] + (-1.0) adp_c::[ADP C10H12N5O10P2] ==> (1.0) pyr_c::[Pyruvate] + (1.0) atp_c::[ATP C10H12N5O13P3]
    raw_map["Supp1"]["PYK"] = ["pykA", "pykF"]
    raw_map["Supp3"]["PYK"] = [
        Dict{String, Any}(
            "raw" => [2970.2, 2833.4, 2641.9, 2780.1, 3020.7, 3346.1, 2889.9, 2283.7, 3005.0, 2180.9, 1678.5, 1517.9], 
            "norm" => [0.231301655, 0.163275862, 0.062317436, 0.135878354, 0.255624487, 0.403222138, 0.191761148, -0.147885281, 0.248106567, -0.214334804, -0.592085887, -0.737181675], 
            "desc" => "MG1655_pykA_b1854  DB_XREF=PID:g1788160  SEG=NC_000913:+1935673,1937115  LEN=1442  DEF=pyruvate kinase II, glucose stimulated",
        ),
        Dict{String, Any}(
            "raw" => [5339.0, 5696.9, 6283.2, 6469.5, 5506.1, 6472.1, 8041.2, 8476.7, 5175.4, 4391.4, 10749.3, 11843.5], 
            "norm" => [-0.334815811, -0.241208278, -0.099885857, -0.057731143, -0.290354548, -0.057151461, 0.256025453, 0.332117369, -0.379714989, -0.616704409, 0.674785449, 0.814638225], 
            "desc" => "MG1655_pykF_b1676  DB_XREF=PID:g1787965  SEG=NC_000913:+1753722,1755134  LEN=1412  DEF=pyruvate kinase I (formerly F), fructose stimulated",
        ),
    ]

    # [95] PYRt2::[Pyruvate transport in via proton symport] Pyruvate transport in via proton symport
    # (-1.0) h_e::[H+] + (-1.0) pyr_e::[Pyruvate] <==> (1.0) h_c::[H+] + (1.0) pyr_c::[Pyruvate]
    raw_map["Supp1"]["PYRt2"] = [] # Not found
    raw_map["Supp3"]["PYRt2"] = [] # Not found

    # [96] RPE::[Ribulose 5-phosphate 3-epimerase] Ribulose 5-phosphate 3-epimerase
    # (-1.0) ru5p__D_c::[D-Ribulose 5-phosphate] <==> (1.0) xu5p__D_c::[D-Xylulose 5-phosphate]
    raw_map["Supp1"]["RPE"] = ["rpe"]
    raw_map["Supp3"]["RPE"] = [
        Dict{String, Any}(
            "raw" => [2872.0, 2703.6, 2764.7, 3043.7, 2506.4, 2639.7, 2550.8, 1900.8, 2074.1, 1853.8, 747.4, 484.8], 
            "norm" => [0.562737508, 0.475563479, 0.5078047, 0.646507927, 0.366298434, 0.441055737, 0.391631546, -0.032711499, 0.093167213, -0.068832636, -1.379365772, -2.003856637], 
            "desc" => "MG1655_rpe_b3386  DB_XREF=PID:g1789788  SEG=NC_000913:-3512020,3512697  LEN=677  DEF=D-ribulose-5-phosphate 3-epimerase",
        )
    ]

    # [97] RPI::[Ribose-5-phosphate isomerase] Ribose-5-phosphate isomerase
    # (-1.0) r5p_c::[Alpha-D-Ribose 5-phosphate] <==> (1.0) ru5p__D_c::[D-Ribulose 5-phosphate]
    raw_map["Supp1"]["RPI"] = ["rpiA", "rpiB"]
    raw_map["Supp3"]["RPI"] = [
        Dict{String, Any}(
            "raw" => [3636.9, 3468.9, 3358.1, 3770.9, 3271.6, 3062.1, 3239.7, 1758.8, 2876.7, 1487.1, 1246.5, 840.7], 
            "norm" => [0.586111481, 0.517880475, 0.471047415, 0.638311115, 0.433398592, 0.337923621, 0.419262447, -0.462006339, 0.247816999, -0.704096112, -0.958714894, -1.526934799], 
            "desc" => "MG1655_rpiA_b2914  DB_XREF=PID:g1789280  SEG=NC_000913:-3056686,3057345  LEN=659  DEF=ribosephosphate isomerase, constitutive",
        ),
        Dict{String, Any}(
            "raw" => [238.2, 209.4, 306.7, 407.1, 399.4, 578.6, 407.1, 996.6, 782.5, 1523.6, 403.0, 397.4], 
            "norm" => [-0.975519874, -1.161431845, -0.610865121, -0.202310066, -0.229858954, 0.304873036, -0.202310066, 1.089321287, 0.740397465, 1.701719001, -0.216913448, -0.237101415], 
            "desc" => "MG1655_rpiB_b4090  DB_XREF=PID:g1790528  SEG=NC_000913:+4310929,4311378  LEN=449  DEF=ribose 5-phosphate isomerase B",
        ),
    ]

    # [98] SUCCt2_2::[Succinate transport via proton symport (2 H)] Succinate transport via proton symport (2 H)
    # (-2.0) h_e::[H+] + (-1.0) succ_e::[Succinate] ==> (2.0) h_c::[H+] + (1.0) succ_c::[Succinate]
    raw_map["Supp1"]["SUCCt2_2"] = [] # Not found
    raw_map["Supp3"]["SUCCt2_2"] = [] # Not found

    # [99] SUCCt3::[Succinate transport out via proton antiport] Succinate transport out via proton antiport
    # (-1.0) h_e::[H+] + (-1.0) succ_c::[Succinate] ==> (1.0) h_c::[H+] + (1.0) succ_e::[Succinate]
    raw_map["Supp1"]["SUCCt3"] = [] # Not found
    raw_map["Supp3"]["SUCCt3"] = [] # Not found

    # [100] SUCDi::[Succinate dehydrogenase (irreversible)] Succinate dehydrogenase (irreversible)
    # (-1.0) q8_c::[Ubiquinone-8] + (-1.0) succ_c::[Succinate] ==> (1.0) q8h2_c::[Ubiquinol-8] + (1.0) fum_c::[Fumarate]
    raw_map["Supp1"]["SUCDi"] = ["sdhABCD"]
    raw_map["Supp3"]["SUCDi"] = [
        Dict{String, Any}(
            "raw" => [6119.3, 5524.7, 5412.1, 7051.7, 7308.7, 7117.2, 5423.7, 1697.9, 2966.0, 1979.5, 420.2, 370.0], 
            "norm" => [1.021218652, 0.873748149, 0.844040519, 1.225823122, 1.277466839, 1.2391618, 0.847129404, -0.828396485, -0.023629379, -0.60701191, -2.84299991, -3.026550801], 
            "desc" => "MG1655_sdhA_b0723  DB_XREF=PID:g1786942  SEG=NC_000913:+755130,756896  LEN=1766  DEF=succinate dehydrogenase, flavoprotein subunit",
        ),
        Dict{String, Any}(
            "raw" => [43.2, 40.7, 44.2, 26.8, 42.6, 54.8, 39.1, 26.9, 35.4, 24.0, 50.6, 42.6], 
            "norm" => [0.181121466, 0.095118948, 0.214136523, -0.507676846, 0.160943584, 0.524266046, 0.037258761, -0.502303674, -0.106160487, -0.666875441, 0.409227538, 0.160943584], 
            "desc" => "CFT073_sdhA_c0801  DB_XREF=26246691  SEG=NC_004431:+782220,783998  LEN=1778  DEF=Succinate dehydrogenase flavoprotein subunit",
        ),
        Dict{String, Any}(
            "raw" => [4457.2, 4375.7, 4514.1, 5328.5, 5561.5, 5286.4, 3928.8, 1443.4, 1749.0, 1225.4, 469.6, 314.6], 
            "norm" => [0.939802962, 0.913179094, 0.958103643, 1.197394728, 1.259139311, 1.185950856, 0.757753992, -0.686863576, -0.409804447, -0.92308198, -2.306830422, -2.88474416], 
            "desc" => "MG1655_sdhB_b0724  DB_XREF=PID:g1786943  SEG=NC_000913:+756912,757628  LEN=716  DEF=succinate dehydrogenase, iron sulfur protein",
        ),
        Dict{String, Any}(
            "raw" => [4510.9, 4485.0, 5305.1, 6469.9, 8035.5, 7105.3, 4667.5, 1804.3, 3275.8, 3221.5, 815.8, 447.8], 
            "norm" => [0.482097555, 0.473790236, 0.716062197, 1.002425665, 1.315070048, 1.137577815, 0.531332274, -0.839878514, 0.020529529, -0.003585154, -1.985030337, -2.850391315], 
            "desc" => "MG1655_sdhC_b0721  DB_XREF=PID:g1786940  SEG=NC_000913:+754400,754789  LEN=389  DEF=succinate dehydrogenase, cytochrome b556",
        ),
        Dict{String, Any}(
            "raw" => [5636.1, 5546.5, 5905.9, 7041.8, 8105.4, 7183.2, 5119.9, 1898.5, 2646.1, 2662.9, 573.0, 418.1], 
            "norm" => [0.807675993, 0.78455646, 0.875135712, 1.128923036, 1.331862168, 1.157605468, 0.669094416, -0.762161217, -0.283153632, -0.274022964, -2.490414172, -2.945101268], 
            "desc" => "MG1655_sdhD_b0722  DB_XREF=PID:g1786941  SEG=NC_000913:+754783,755130  LEN=347  DEF=succinate dehydrogenase, hydrophobic subunit",
        ),
    ]

    # [101] SUCOAS::[Succinyl-CoA synthetase (ADP-forming)] Succinyl-CoA synthetase (ADP-forming)
    # (-1.0) succ_c::[Succinate] + (-1.0) atp_c::[ATP C10H12N5O13P3] + (-1.0) coa_c::[Coenzyme A] <==> (1.0) pi_c::[Phosphate] + (1.0) succoa_c::[Succinyl-CoA] + (1.0) adp_c::[ADP C10H12N5O10P2]
    raw_map["Supp1"]["SUCOAS"] = ["sucCD"]
    raw_map["Supp3"]["SUCOAS"] = [
        Dict{String, Any}(
            "raw" => [8574.7, 8609.9, 8198.4, 8591.5, 9305.5, 10898.2, 9953.1, 5342.3, 6941.5, 4505.6, 1804.1, 646.1], 
            "norm" => [0.605992944, 0.611903228, 0.541249127, 0.608816782, 0.723990417, 0.951924713, 0.821052685, -0.07663226, 0.301154197, -0.322374014, -1.642813945, -3.124263873], 
            "desc" => "MG1655_sucC_b0728  DB_XREF=PID:g1786948  SEG=NC_000913:+762237,763403  LEN=1166  DEF=succinyl-CoA synthetase, beta subunit",
        ),
        Dict{String, Any}(
            "raw" => [7701.5, 7844.2, 7124.8, 8609.7, 9392.6, 9975.2, 10869.2, 6345.5, 7186.0, 5131.3, 1760.2, 1010.5], 
            "norm" => [0.399416271, 0.425903128, 0.287126325, 0.560239777, 0.685801379, 0.772622576, 0.896450662, 0.120010655, 0.299465744, -0.186398817, -1.72998383, -2.53065387], 
            "desc" => "MG1655_sucD_b0729  DB_XREF=PID:g1786949  SEG=NC_000913:+763403,764272  LEN=869  DEF=succinyl-CoA synthetase, alpha subunit",
        ),
    ]

    # [102] TALA::[Transaldolase] Transaldolase
    # (-1.0) s7p_c::[Sedoheptulose 7-phosphate] + (-1.0) g3p_c::[Glyceraldehyde 3-phosphate] <==> (1.0) e4p_c::[D-Erythrose 4-phosphate] + (1.0) f6p_c::[D-Fructose 6-phosphate]
    raw_map["Supp1"]["TALA"] = ["talB"]
    raw_map["Supp3"]["TALA"] = [
        Dict{String, Any}(
            "raw" => [9249.2, 10278.6, 9947.8, 12955.0, 10029.1, 9369.7, 10395.9, 4039.5, 7129.4, 5906.6, 2779.6, 1807.3], 
            "norm" => [0.43834696, 0.590590243, 0.543395875, 0.924455483, 0.555138614, 0.457021229, 0.606961129, -0.756804897, 0.06279904, -0.208653712, -1.296104341, -1.917145623], 
            "desc" => "MG1655_talB_b0008  DB_XREF=PID:g1786189  SEG=NC_000913:+8238,9191  LEN=953  DEF=transaldolase B",
        ),
    ]

    # [103] THD2::[NAD(P) transhydrogenase] NAD(P) transhydrogenase
    # (-2.0) h_e::[H+] + (-1.0) nadh_c::[Nicotinamide adenine dinucleotide - reduced] + (-1.0) nadp_c::[Nicotinamide adenine dinucleotide phosphate] ==> (2.0) h_c::[H+] + (1.0) nad_c::[Nicotinamide adenine dinucleotide] + (1.0) nadph_c::[Nicotinamide adenine dinucleotide phosphate - reduced]
    raw_map["Supp1"]["THD2"] = ["pntAB"]
    raw_map["Supp3"]["THD2"] = [
        Dict{String, Any}(
            "raw" => [13071.5, 14675.0, 14760.1, 18214.6, 10825.3, 10443.5, 11896.1, 5841.4, 9354.4, 8085.3, 3590.7, 2353.8], 
            "norm" => [0.548844622, 0.715780421, 0.724122413, 1.027515231, 0.276826924, 0.22502521, 0.412908597, -0.613193998, 0.066136943, -0.144206872, -1.315243055, -1.924516436], 
            "desc" => "MG1655_pntA_b1603  DB_XREF=PID:g1787887  SEG=NC_000913:-1674395,1675927  LEN=1532  DEF=pyridine nucleotide transhydrogenase, alpha subunit",
        ),
        Dict{String, Any}(
            "raw" => [9321.5, 9529.0, 10441.3, 10384.0, 6494.6, 6282.7, 6642.8, 1982.1, 4586.7, 3669.2, 1497.2, 993.0], 
            "norm" => [0.974903126, 1.006665818, 1.138570438, 1.13063138, 0.45358167, 0.405725689, 0.486132475, -1.258629253, -0.048202455, -0.370193459, -1.66339205, -2.255793381], 
            "desc" => "MG1655_pntB_b1602  DB_XREF=PID:g1787886  SEG=NC_000913:-1672996,1674384  LEN=1388  DEF=pyridine nucleotide transhydrogenase, beta subunit",
        ),
    ]

    # [104] TKT1::[Transketolase] Transketolase
    # (-1.0) r5p_c::[Alpha-D-Ribose 5-phosphate] + (-1.0) xu5p__D_c::[D-Xylulose 5-phosphate] <==> (1.0) s7p_c::[Sedoheptulose 7-phosphate] + (1.0) g3p_c::[Glyceraldehyde 3-phosphate]
    raw_map["Supp1"]["TKT1"] = ["tktA"]
    raw_map["Supp3"]["TKT1"] = [
        Dict{String, Any}(
            "raw" => [6794.5, 6776.0, 6623.7, 6861.5, 7052.1, 4380.7, 5070.7, 1113.4, 3926.0, 1985.5, 610.1, 475.7], 
            "norm" => [1.126439225, 1.122505712, 1.089709169, 1.140595837, 1.180124769, 0.493233256, 0.704256759, -1.482956175, 0.33513201, -0.648425802, -2.350810527, -2.709804233], 
            "desc" => "MG1655_tktA_b2935  DB_XREF=PID:g2367177  SEG=NC_000913:-3077663,3079654  LEN=1991  DEF=transketolase 1 isozyme",
        ),
    ]

    # [105] TKT2::[Transketolase] Transketolase
    # (-1.0) xu5p__D_c::[D-Xylulose 5-phosphate] + (-1.0) e4p_c::[D-Erythrose 4-phosphate] <==> (1.0) f6p_c::[D-Fructose 6-phosphate] + (1.0) g3p_c::[Glyceraldehyde 3-phosphate]
    raw_map["Supp1"]["TKT2"] = ["tktB"]
    raw_map["Supp3"]["TKT2"] = [
        Dict{String, Any}(
            "raw" => [2263.4, 2325.7, 2460.8, 2527.2, 1343.0, 3012.0, 2252.6, 4102.8, 2687.9, 1617.8, 5242.7, 5742.6], 
            "norm" => [-0.262319406, -0.223145963, -0.141683565, -0.103271132, -1.015351669, 0.149910796, -0.269219821, 0.595797855, -0.014331508, -0.746777708, 0.94949902, 1.080893101], 
            "desc" => "MG1655_tktB_b2465  DB_XREF=PID:g1788808  SEG=NC_000913:+2577656,2579659  LEN=2003  DEF=transketolase 2 isozyme",
        ),
        Dict{String, Any}(
            "raw" => [117.4, 113.3, 138.5, 111.5, 111.4, 131.9, 149.5, 179.0, 123.3, 176.8, 204.7, 328.0], 
            "norm" => [-0.34234316, -0.393627707, -0.103889592, -0.416731858, -0.418026335, -0.174331003, 0.006369916, 0.266184019, -0.271602768, 0.248342707, 0.459735534, 1.139920247], 
            "desc" => "CFT073_tktB_c2990  DB_XREF=26248832  SEG=NC_004431:+2850331,2852397  LEN=2066  DEF=Transketolase 2",
        )
    ]

    # [106] TPI::[Triose-phosphate isomerase] Triose-phosphate isomerase
    # (-1.0) dhap_c::[Dihydroxyacetone phosphate] <==> (1.0) g3p_c::[Glyceraldehyde 3-phosphate]
    raw_map["Supp1"]["TPI"] = ["tpiA"]
    raw_map["Supp3"]["TPI"] = [
        Dict{String, Any}(
            "raw" => [7587.7, 9408.8, 9331.9, 10843.5, 7596.8, 8720.2, 10172.3, 10094.8, 8915.7, 6217.5, 9980.4, 8687.6], 
            "norm" => [-0.224771386, 0.085576708, 0.073736823, 0.290324566, -0.223042185, -0.024072801, 0.198139985, 0.187106398, 0.007914048, -0.512099423, 0.170663613, -0.029476346], 
            "desc" => "MG1655_tpiA_b3919  DB_XREF=PID:g1790353  SEG=NC_000913:-4108320,4109087  LEN=767  DEF=triosephosphate isomerase",
        ),
    ]

    # [107] UDPG4E::[UDPglucose 4-epimerase] Alternate Carbon Metabolism
    # (-1.0) udpg_c::[UDPglucose] <==> (1.0) udpgal_c::[UDPgalactose]
    raw_map["Supp1"]["UDPG4E"] = ["galE"]
    raw_map["Supp3"]["UDPG4E"] = [
        Dict{String, Any}(
            "raw" => [3248.8, 3439.4, 6224.6, 11543.8, 11162.7, 8508.8, 6151.1, 2194.7, 991.6, 1549.3, 1133.0, 836.1], 
            "norm" => [-0.00966567, 0.072584307, 0.928408528, 1.819473702, 1.771041516, 1.37938308, 0.911271828, -0.575548856, -1.721742426, -1.077956074, -1.529424741, -1.967825194], 
            "desc" => "MG1655_galE_b0759  DB_XREF=PID:g1786974  SEG=NC_000913:-790262,791278  LEN=1016  DEF=UDP-galactose-4-epimerase",
        )
    ]

    # [108] UGLT::[UDPglucose--hexose-1-phosphate uridylyltransferase] Alternate Carbon Metabolism
    # (-1.0) gal1p_c::[Alpha-D-Galactose 1-phosphate] + (-1.0) udpg_c::[UDPglucose] <==> (1.0) udpgal_c::[UDPgalactose] + (1.0) g1p_c::[D-Glucose 1-phosphate]
    raw_map["Supp1"]["UGLT"] = ["galU"]
    raw_map["Supp3"]["UGLT"] = [
        Dict{String, Any}(
            "raw" => [1689.2, 2259.4, 2506.0, 2167.4, 1786.9, 2237.1, 2946.9, 2974.5, 2974.4, 1392.5, 2373.7, 2777.3], 
            "norm" => [-0.43329837, -0.013698817, 0.135747892, -0.073673091, -0.352179623, -0.028008775, 0.369559582, 0.383008656, 0.382960153, -0.711961195, 0.057499089, 0.284044501], 
            "desc" => "MG1655_galU_b1236  DB_XREF=PID:g1787488  SEG=NC_000913:+1290680,1291588  LEN=908  DEF=glucose-1-phosphate uridylyltransferase",
        ),
    ]

    ## ------------------------------------------------
    # build a matriz
    @stage! exp_mat = Dict()
    for (rxn, genes) in raw_map["Supp3"]
        isempty(genes) && continue
        nrows = length(genes)
        ncols = length(genes[1]["raw"])
        exp_mat[rxn] = zeros(nrows, ncols)
        for (ri, obj) in enumerate(genes)
            exp_mat[ri, :] .= obj["raw"]
        end
    end

    exp_mat
end

## ------------------------------------------------
## ------------------------------------------------
## ------------------------------------------------
let
    f = Figure()
    ax = Axis(f[1,1]; xlabel = "time", ylabel = "gen expression")
    rxns = ["ACALD", "PTAr", "UDPG4E", "ALCD2x", "PDH", "CO2t", "PYK", "EX_nh4_e", "MALt2_2", "CS", "GLYCtpp", "G3PD2", "PGM", "TKT1", "EX_mal__L_e", "ACONTa", "EX_pi_e", "GLNS", "ICL", "PGMT", "EX_o2_e", "FBA", "EX_gln__L_e", "EX_glc__D_e", "FORt2", "SUCCt3", "G6PDH2r", "EX_malt_e", "AKGDH", "TKT2", "FRD7", "SUCOAS", "BIOMASS_Ecoli_core_w_GAM", "FBP", "ICDHyr", "AKGt2r", "GLUSy", "TPI", "GALKr", "FORt", "ACONTb", "EX_ac_e", "GLNabc", "EX_akg_e", "EX_fru_e", "RPE", "ACKr", "EX_glyc_e", "THD2", "D_LACt2", "EX_glu__L_e", "PFL", "RPI", "TALA", "ATPM", "ACt2r", "EX_etoh_e", "PPCK", "NH4t", "PGL", "NADTRHD", "PGK", "GALabcpp", "LDH_D", "ME1", "PIt2r", "EX_h2o_e", "EX_succ_e", "ATPS4r", "EX_acald_e", "PYRt2", "EX_h_e", "UGLT", "GLCpts", "GLUDy", "CYTBD", "EX_gal_e", "FUMt2_2", "FRUpts2", "GAPD", "H2Ot", "NADH16", "PPC", "GLYK", "EX_for_e", "PFK", "MDH", "PGI", "O2t", "ME2", "EX_pyr_e", "EX_co2_e", "AMALT1", "GND", "GLUN", "SUCCt2_2", "EX_fum_e", "ETOHt2r", "MALTabcpp", "ADK1", "ACALDt", "EX_lac__D_e", "SUCDi", "ENO", "MALS", "GLUt2r", "PPS", "FUM"]
    for (rxn, genes) in raw_map["Supp3"]
        # rxn == rand(rxns) || continue
        rxn == rand(rxns) || continue
        @show rxn
        exp_sum = nothing
        genid = ""
        for gene in genes
            genid = gene["desc"][1:20]
            exp_sum = isnothing(exp_sum) ? gene["raw"] : exp_sum + gene["raw"]
            # lines!(ax, patt; label = genid)
            # lines!(ax, gene["raw"]; label = genid)
        end
        exp_patt = exp_sum .- minimum(exp_sum)
        exp_patt = exp_patt ./ maximum(exp_patt)
        lines!(ax, exp_patt; label = genid)
    end
    axislegend(ax)
    f
end