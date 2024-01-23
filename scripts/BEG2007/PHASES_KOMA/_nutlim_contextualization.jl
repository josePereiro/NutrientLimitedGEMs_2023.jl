@time let 
    using MetX
    using Clp
    using DataFrames
end

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
include("1_setup.jl")
include("2_utils.jl")

## --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 

# load network
net0 = pull_net("iJO1366")
size(net0) # (1805, 2583)

## --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
let
    for (modelid, net0) = [
            ("core", pull_net("ecoli_core_Beg2007")), 
            ("iJO1366", pull_net("iJO1366"))
        ]
        global results = DataFrame()
        
        # contextualize (from experimental condition)
        lb!(net0, "EX_glc__D_e", -5.0) # [mmol/ gDW h]
        lb!(net0, "EX_lac__D_e", -5.0) # [mmol/ gDW h]
        lb!(net0, "EX_malt_e", -5.0)   # [mmol/ gDW h]
        lb!(net0, "EX_gal_e", -5.0)    # [mmol/ gDW h]
        lb!(net0, "EX_glyc_e", -5.0)   # [mmol/ gDW h]
        lb!(net0, "EX_ac_e", 0.0)      # [mmol/ gDW h]
        ub!(net0, "EX_ac_e", 5.0)      # [mmol/ gDW h]

        # 90% biom/nut shadow prices
        opm = FBAOpModel(net0, Clp.Optimizer)
        set_linear_obj!(opm, extras(net0, "BIOM"), MAX_SENSE)
        for exch in [
                "EX_glc__D_e", "EX_lac__D_e", "EX_malt_e",
                "EX_gal_e", "EX_glyc_e"
            ]
            exch in colids(opm) || continue
            test_points = lb(opm, exch) .* [1.0, 0.95, 0.9]
            price, price_err, _, _ = bound_dual_prices(
                opm, exch, test_points, :lb; 
                dovars = false
            )
            push!(results, (;exch, price, price_err))
        end 

        println("\n", "-"^50)
        println(modelid, " ", size(net0))
        println("biomass ", objective_value(optimize!(opm)))
        println(results)
    end
end 
## --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
