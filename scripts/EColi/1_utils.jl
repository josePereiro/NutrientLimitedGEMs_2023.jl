## ------------------------------------------------------------------
function _setup_ecoli(;solver = LP_SOLVER)

    net = pull_net("ecoli_core")
    clampbounds!(net, -1000.0, 1000.0)

    # q_nh4 = 3.24 # (mmol/ gCDW h)
    q_nh4 = 4.0 # (mmol/ gCDW h)
    # q_glc = 10.5 # (mmol/ gCDW h)
    q_glc = 10.0 # (mmol/ gCDW h)
    @show q_nh4
    @show q_glc

    ex_glc = extras(net, "EX_GLC")
    ex_nh4 = extras(net, "EX_NH4")
    lb!(net, ex_nh4, -q_nh4)
    lb!(net, ex_glc, -q_glc)

    # max biom
    opm = FBAFluxOpModel(net, solver)
    optimize!(opm)
    biom_id = extras(net, "BIOM")
    biom0 = solution(opm, biom_id)
    @show biom0
    
    return box(net, solver)

end

## ------------------------------------------------------------------
return nothing