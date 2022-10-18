@time begin
    using COBREXA
    # using NutrientLimitedGEMs
    # const NL = NutrientLimitedGEMs
    # using GLPK
    # using ProjAssistant
    # using Plots
    # using Plots.Measures
    # using JuMP
    # using Ipopt
    using MetNets
end

## ------------------------------------------------------------------
let
    # global net = NL.load_toy_model4()
    # target = "Ex_s"
    # biom = "biomass"
    
    global net = load_ecoli()
    target = "R_EX_gln__L_e"
    biom = "R_BIOMASS_Ecoli_core_w_GAM"
    
    biomidx = rxnidx(net, biom)
    change_objective!(net, biom)

    bounds!(net, target; lb = -1.0)
    bounds!(net, biom; ub = 0.5)
    global dense0 = dense_fba(net)
    global fba0 = flux_balance_analysis_vec(net, Ipopt.Optimizer)
    
    bounds!(net, target; lb = -0.5)
    global dense1 = dense_fba(net)
    global fba1 = flux_balance_analysis_vec(net, Ipopt.Optimizer)

    println("\n\n\n", "-"^60)
    println(rpad("id", 25), "\tf0\tf1\tdiff\td0\td1\tdiff", "\n")
    digits = 2
    for (id, d0, d1, f0, f1) in zip(net.rxns, dense0, dense1, fba0, fba1)
        println(
            rpad(id, 25), "\t", 
            round(f0; digits), "\t", 
            round(f1; digits), "\t", 
            round(f0 - f1; digits), "\t",
            round(d0; digits), "\t", 
            round(d1; digits), "\t", 
            round(d0 - d1; digits), "\t", 
        )
    end
    println(
        "Dense nz count ", 
        count(.!isapprox.(dense0 .- dense1, 0.0; atol = 1e-2))
    )
    println(
        "Fba nz count ", 
        count(.!isapprox.(fba1 .- fba0, 0.0; atol = 1e-2))
    )
end
