@time begin
    using MetXBase
    using MetXGEMs
    using MetXOptim
    using MetXNetHub
    using NutrientLimitedGEMs
    using Clp
    using LinearAlgebra
    using CairoMakie
    # using COBREXA
end

## ------------------------------------------------------------
# Compute ctx shadow price
let
    net = pull_net("ecoli_core_Beg2007")
    biom_id = extras(net, "BIOM")
    linear_weights!(net, biom_id, 1.0)

    exchs = [
        "EX_glc__D_e", "EX_lac__D_e", "EX_malt_e",
        "EX_gal_e", "EX_glyc_e", "EX_ac_e", "EX_ac_e"
    ]
    
    f = Figure()
    ax = Axis(f[1,1]; 
        xlabel = "var index (sorted)", 
        ylabel = "lb shadow prices"
    )
    for exch in exchs
        opm = FBAOpModel(net, Clp.Optimizer)
        test_points = lb(net, exch) .* [1.0, 0.95, 0.9]
        obj_m, obj_err, vars_ms, vars_errs = bound_dual_prices(
            opm, exch, test_points, :lb; 
            dovars = true
        )
        sperm = sortperm(vars_ms)
        xs = eachindex(sperm)
        ys = vars_ms[sperm]
        errs = vars_errs[sperm]
        lines!(ax, xs, ys; label = exch, linewidth = 4)
        errorbars!(ax, xs, ys, errs)
    end
    axislegend(ax; position = :rb)
    f
end

## ------------------------------------------------------------