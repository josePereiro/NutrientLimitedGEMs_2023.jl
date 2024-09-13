## ------------------------------------------------------------------
function _setup_iJR904(;solver = LP_SOLVER)

    net = pull_net("iJR904")

    # close all intakes, open outtakes
    for rxn in net.rxns
        startswith(rxn, "R_EX_") || continue
        bounds!(net, rxn, 0.0, 1000.0)
    end

    # Medium
    q_glc = 10.5 # (mmol/ gCDW h)
    q_nh4 = 3.24 # (mmol/ gCDW h)
    
    lb!(net, "R_EX_glc__D_e", -q_glc)
    lb!(net, "R_EX_nh4_e", -q_nh4)
    lb!(net, "R_EX_fe2_e", -1000.0)
    lb!(net, "R_EX_h2o_e", -1000.0)
    lb!(net, "R_EX_h_e", -1000.0)
    lb!(net, "R_EX_k_e", -1000.0)
    lb!(net, "R_EX_na1_e", -1000.0)
    lb!(net, "R_EX_o2_e", -1000.0)
    lb!(net, "R_EX_pi_e", -1000.0)
    lb!(net, "R_EX_so4_e", -1000.0)
    
    # # fva
    # cid = (hash(net), hash(solver))
    # net = lcache(NutrientLimitedGEMs_2023, cid) do 
    #    box(net, solver; verbose = true)
    # end

    # max biom
    opm = FBAOpModel(net, solver)
    optimize!(opm)
    biom_id = extras(net, "BIOM")
    biom0 = solution(opm, biom_id)
    @show biom0

    return net

end

## ------------------------------------------------------------------
return nothing