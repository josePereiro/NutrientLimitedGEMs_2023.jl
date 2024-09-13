## ------------------------------------------------------------------
@time begin
    using NutrientLimitedGEMs_2023
    const NL = NutrientLimitedGEMs_2023
    
    using ProjFlows
    using ContextDBs
    using MetX
    using Gurobi
end

## ------------------------------------------------------------------
let
    
end
## ------------------------------------------------------------------
## ------------------------------------------------------------------
## ------------------------------------------------------------------
# Experimental data
_, (CalDat, _) = lprocdat(PROJ,
    ["Calzadilla_et_al"], "Calzadilla.bundle", ".jls"; 
    verbose = true
)
nothing


# ------------------------------------------------------------------
let
    st = 5
    solver = Gurobi.Optimizer

    model_id = "Martinez_Monge_HEK293"
    net = pull_net(model_id)
    clampbounds!(net, -1000.0, 1000.0)
    
    # contextualization
    # lb = - c D / Xv

    stdat = CalDat[st]
    D = stdat["D"]["val"]
    Xv = stdat["Xv"]["val"]

    # ["ALA", "ARG", "ASN", "ASP", "GLN", "GLU", "GLC", 
    # "GLY", "ILE", "LAC", "LEU", "LYS", "MET", "PHE", 
    # "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
    for (CalId, HekID) in [
        ("ALA", "R_EX_ala_L_LPAREN_e_RPAREN_"),
        # ("ARG", "R_EX_arg_L_LPAREN_e_RPAREN_"), # TODO: ask Calzadilla
        ("ASN", "R_EX_asn_L_LPAREN_e_RPAREN_"),
        ("ASP", "R_EX_asp_L_LPAREN_e_RPAREN_"),
        ("GLN", "R_EX_gln_L_LPAREN_e_RPAREN_"),
        ("GLU", "R_EX_glu_L_LPAREN_e_RPAREN_"),
        ("GLC", "R_EX_glc_LPAREN_e_RPAREN_"),
        ("GLY", "R_EX_gly_LPAREN_e_RPAREN_"),
        ("ILE", "R_EX_ile_L_LPAREN_e_RPAREN_"),
        ("LAC", "R_EX_lac_L_LPAREN_e_RPAREN_"),
        ("LEU", "R_EX_leu_L_LPAREN_e_RPAREN_"),
        ("LYS", "R_EX_lys_L_LPAREN_e_RPAREN_"),
        ("MET", "R_EX_met_L_LPAREN_e_RPAREN_"),
        ("PHE", "R_EX_phe_L_LPAREN_e_RPAREN_"),
        ("PRO", "R_EX_pro_L_LPAREN_e_RPAREN_"),
        ("SER", "R_EX_ser_L_LPAREN_e_RPAREN_"),
        ("THR", "R_EX_thr_L_LPAREN_e_RPAREN_"),
        ("TRP", "R_EX_trp_L_LPAREN_e_RPAREN_"),
        ("TYR", "R_EX_tyr_L_LPAREN_e_RPAREN_"),
        ("VAL", "R_EX_val_L_LPAREN_e_RPAREN_"),
    ]
        c = stdat["c$(CalId)"]["val"]
        lb = -(c * D / Xv)
        
        println(CalId, " lb: ", lb)
        lb!(net, HekID, lb)

    end
    
    # TODO: ask Calzadilla
    lb!(net, "R_EX_his_L_LPAREN_e_RPAREN_", -1000.0)
    lb!(net, "R_EX_nh4_LPAREN_e_RPAREN_", -1000.0)
    
    lb!(net, "R_EX_o2_LPAREN_e_RPAREN_", -1000.0)
    lb!(net, "R_EX_h_LPAREN_e_RPAREN_", -1000.0)
    lb!(net, "R_EX_pi_LPAREN_e_RPAREN_", -1000.0)
    lb!(net, "R_EX_h2o_LPAREN_e_RPAREN_", -1000.0)
    lb!(net, "R_EX_co2_LPAREN_e_RPAREN_", -1000.0)
    lb!(net, "R_EX_hco3_LPAREN_e_RPAREN_", -1000.0)

    # max biom
    opm = FBAOpModel(net, solver)
    optimize!(opm)
    biom_id = extras(net, "BIOM")
    biom0 = solution(opm, biom_id)
    @show biom0
    @show D

    nothing
end