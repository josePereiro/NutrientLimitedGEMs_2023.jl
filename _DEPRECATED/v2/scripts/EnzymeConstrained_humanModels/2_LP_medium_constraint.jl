## ------------------------------------------------------------------
@time begin
    using MetXBase
    using MetXOptim
    using MetXNetHub
    using MetXCultureHub
    using Gurobi
    using ContextDBs
    using ExtractMacro
    using Plots
    using Serialization
end

# ------------------------------------------------------------------
include("1_utils.jl")

# ------------------------------------------------------------------
let
    global cul_db = pull_cul("rathCharacterisationCellGrowth2017");
    global cul_iters = query(cul_db, ["ROOT", "Iters"])

    gemids = [
        # "HOP62", 
        # "HS_578T",
        # "MDMAMB_231", 
        # "NCI_H226",  
        "UO_31"
    ]
    # gemids = ["HOP62", "HOP92"]
    culids = ["A", "B", "C", "D", "E"]
    # biom_mod = "AGE1.HN.AAT" # Rath study the producer line
    biom_mod = nothing # Rath study the producer line

    p = plot(; xlabel = "exp growth [1/h]", ylabel = "LP max growth [1/h]")

    for gemid in gemids
        
        println("\n", "="^60)
        @show gemid
        global net = pull_net("SysBioChalmers_EnzymeConstrained_humanModels", gemid, biom_mod)
        @show size(net)
        _test_fba(net)
        println()
        
        exp_bioms = []
        max_bioms = []
        for culid in culids
            println("\n", "."^60)
            @show culid
            _set_rath_medium_bounds!(net, cul_db, culid)

            # culture
            stst_db = queryall(cul_db, ["ROOT", "culid" => culid])
            exp_biom = query(stst_db, ["ROOT", "apiid" => "μ"])["val"]

            # LP
            opm = FBAOpModel(net, Gurobi.Optimizer)
            optimize!(opm)
            biom =  solution(opm, "biomass_human")
            @show biom
            @show exp_biom
            push!(exp_bioms, exp_biom)
            push!(max_bioms, biom)
        end

        scatter!(p, exp_bioms, max_bioms; label = gemid)
    end

    # plot!(; 
    #     xlim = [-0.0, 0.018], 
    #     ylim = [-0.0, 0.018]
    # )
    p
end

## ------------------------------------------------------------------
let
    # culture data
    culid = "C"
    stst_db = queryall(cul_db, ["ROOT", "culid" => culid])
    c_db = queryall(stst_db, ["ROOT", "apiid" => r"c_"])
    @assert all(c_db[:, "unit"] .== "mM")
    D_db = query(stst_db, ["ROOT", "apiid" => "D"])
    @show D = D_db["val"]
    Xv_db = query(stst_db, ["ROOT", "apiid" => "Xv"])
    @show Xv = Xv_db["val"]
    nothing
end