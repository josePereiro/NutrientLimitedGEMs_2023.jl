const ECOLI_BIOM_ID = "R_BIOMASS_Ecoli_core_w_GAM"
const ECOLI_EX_NH4_ID = "R_EX_nh4_e_bkwd"
const ECOLI_EX_GLC_ID = "R_EX_glc__D_e_bkwd"

function load_ecoli(; dobox = false)

    # model = MetNets.ecoli_core_model()
    fn = rawdir(@__MODULE__, "e_coli_core.xml")
    model = MetNets.read_model(fn)

    # pre-process
    # open all exchanges
    MetNets.ub!.([model], MetNets.exchanges(model), 100.0)

    # open gln transportation
    MetNets.ub!(model, "R_GLNabc", 100.0)

    return dobox ? MetLP.fva_preprocess(model; verbose = false) : model
end

function ecoli_core_test_model()
    
    # net
    net = load_ecoli(; dobox = false)
    
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
        MetNets.bounds!(net, rxn; ub = 0.0)
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

    MetNets.bounds!(net, "R_EX_glc__D_e_bkwd"; ub = q_glc)
    MetNets.bounds!(net, "R_EX_nh4_e_bkwd"; ub = q_nh4)
    MetNets.bounds!(net, "R_EX_o2_e_bkwd"; ub = 1000.0)
    MetNets.bounds!(net, "R_EX_h2o_e_bkwd"; ub = 1000.0)
    MetNets.bounds!(net, "R_EX_pi_e_bkwd"; ub = 1000.0)
    
    # under bound biomass
    objid = "R_BIOMASS_Ecoli_core_w_GAM"
    sol = MetLP.fba(net, objid)
    objval = MetNets.av(net, sol, objid)
    @show objval
    MetNets.bounds!(net, objid; ub = objval * 2)
    
    return net

end