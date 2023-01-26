## ------------------------------------------------------------------
# Martinez_Monge_HEK293
function _setup_heknet(;solver = LP_SOLVER)

    model_id = "Martinez_Monge_HEK293"
    net = pull_net(model_id)
    clampbounds!(net, -1000.0, 1000.0)
    
    # contextualization
    lb!(net, "R_EX_tyr_L_LPAREN_e_RPAREN_", -1.0)
    lb!(net, "R_EX_his_L_LPAREN_e_RPAREN_", -1.0)
    lb!(net, "R_EX_asn_L_LPAREN_e_RPAREN_", -1.0)
    lb!(net, "R_EX_arg_L_LPAREN_e_RPAREN_", -1.0)
    lb!(net, "R_EX_ile_L_LPAREN_e_RPAREN_", -1.0)
    lb!(net, "R_EX_phe_L_LPAREN_e_RPAREN_", -1.0)
    lb!(net, "R_EX_met_L_LPAREN_e_RPAREN_", -1.0)
    lb!(net, "R_EX_thr_L_LPAREN_e_RPAREN_", -1.0)
    lb!(net, "R_EX_lys_L_LPAREN_e_RPAREN_", -1.0)
    lb!(net, "R_EX_leu_L_LPAREN_e_RPAREN_", -1.0)
    lb!(net, "R_EX_val_L_LPAREN_e_RPAREN_", -1.0)
    lb!(net, "R_EX_trp_L_LPAREN_e_RPAREN_", -1.0)
    
    lb!(net, "R_EX_glc_LPAREN_e_RPAREN_", -1.0)
    
    lb!(net, "R_EX_glu_L_LPAREN_e_RPAREN_", -1.0)
    lb!(net, "R_EX_gln_L_LPAREN_e_RPAREN_", -1.0)
    # lb!(net, "R_EX_lac_L_LPAREN_e_RPAREN_", -1000.0)
    # lb!(net, "R_EX_nh4_LPAREN_e_RPAREN_", -1000.0)
    
    lb!(net, "R_EX_o2_LPAREN_e_RPAREN_", -1000.0)
    lb!(net, "R_EX_h_LPAREN_e_RPAREN_", -1000.0)
    lb!(net, "R_EX_pi_LPAREN_e_RPAREN_", -1000.0)
    lb!(net, "R_EX_h2o_LPAREN_e_RPAREN_", -1000.0)
    lb!(net, "R_EX_co2_LPAREN_e_RPAREN_", -1000.0)
    lb!(net, "R_EX_hco3_LPAREN_e_RPAREN_", -1000.0)

    # max biom
    opm = FBAFluxOpModel(net, solver)
    optimize!(opm)
    biom_id = extras(net, "BIOM")
    biom0 = solution(opm, biom_id)
    @show biom0

    return net
end

