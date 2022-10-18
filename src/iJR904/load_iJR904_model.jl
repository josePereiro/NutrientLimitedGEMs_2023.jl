const iJR904_BIOM_ID = "BiomassEcoli"
const iJR904_EX_NH4_ID = "EX_nh4_LPAREN_e_RPAREN__bkwd"
const iJR904_EX_GLC_ID = "EX_glc_LPAREN_e_RPAREN__bkwd"

function load_iJR904_model(; clearcache = false)

    # cache
    cfile = procdir(NutrientLimitedGEMs, "iJR904_Folsom2015.jls")
    clearcache && rm(cfile)
    isfile(cfile) && return deserialize(cfile)::MetNets.MetNet

    # net
    mfile = rawdir(NutrientLimitedGEMs, "iJR904.mat")
    @assert isfile(mfile)
    net = MetNets.read_model(mfile)
    
    # open all exchanges
    exchs = MetNets.exchanges(net)
    MetNets.bounds!.([net], exchs; lb = -10000.0, ub = 10000.0)
    
    # foward_defined_model
    net = MetNets.foward_defined_model(net)
    
    # close all intakes
    exchs = MetNets.exchanges(net)
    for rxni in exchs
        rxn = net.rxns[rxni]
        contains(rxn, "EX_") || continue
        contains(rxn, "_bkwd") || continue
        MetNets.bounds!(net, rxni; ub = 0.0)
    end
    
    # open medium
    # medium from 10.1099/mic.0.000118 
    # Supp mat: 000118s3 (biomass)
    # Supp mat: 000118s4 (fluxes, yield, D)

    # TODO: Move to Folsom2015
    # From Y_gX/gN (g/g) * MM_N (g/mol) * c_nh4 (mM) * 1e-3 = X (g/L)
    Xs = [8.81, 8.82, 9.84, 9.22] .* 14 .* 1.31 .* 1e-3

    # From D (1/h) * c_nh4 (mM) / max(X) (gCDW/L) = q_nh4 (mmol/ gCDW h)
    q_nh4 = [0.1, 0.2, 0.3, 0.4] .* 1.31 ./ maximum(Xs)
    q_nh4 = maximum(q_nh4)

    q_glc = 10.5 # (mmol/ gCDW h)

    # Medium
    MetNets.bounds!(net, "EX_glc_LPAREN_e_RPAREN__bkwd"; ub = q_glc)
    MetNets.bounds!(net, "EX_nh4_LPAREN_e_RPAREN__bkwd"; ub = q_nh4)
    MetNets.bounds!(net, "EX_fe2_LPAREN_e_RPAREN__bkwd"; ub = 1000.0)
    MetNets.bounds!(net, "EX_h2o_LPAREN_e_RPAREN__bkwd"; ub = 1000.0)
    MetNets.bounds!(net, "EX_h_LPAREN_e_RPAREN__bkwd"; ub = 1000.0)
    MetNets.bounds!(net, "EX_k_LPAREN_e_RPAREN__bkwd"; ub = 1000.0)
    MetNets.bounds!(net, "EX_na1_LPAREN_e_RPAREN__bkwd"; ub = 1000.0)
    MetNets.bounds!(net, "EX_o2_LPAREN_e_RPAREN__bkwd"; ub = 1000.0)
    MetNets.bounds!(net, "EX_pi_LPAREN_e_RPAREN__bkwd"; ub = 1000.0)
    MetNets.bounds!(net, "EX_so4_LPAREN_e_RPAREN__bkwd"; ub = 1000.0)
    
    # under bound biomass
    objid = "BiomassEcoli"
    sol = MetLP.fba(net, objid)
    objval = MetNets.av(net, sol, objid)
    @show objval
    MetNets.bounds!(net, objid; ub = objval * 2)

    # save
    serialize(cfile, net)
    
    return net

end