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

    exchs = [
        "EX_glc__D_e", "EX_lac__D_e", "EX_malt_e",
        "EX_gal_e", "EX_glyc_e", "EX_ac_e", "EX_ac_e"
    ]

    f = Figure(resolution = (1000, 1000))
    for (i, netid) in enumerate(["ecoli_core_Beg2007", "iJO1366"])

        net = pull_net(netid)
        # net = pull_net("ecoli_core_Beg2007")
        # net = pull_net("iJO1366")
        biom_id = extras(net, "BIOM")
        linear_weights!(net, biom_id, 1.0)

        gf = 0.18 
        lb!(net, "EX_glc__D_e", -0.9 * 60 * gf)
        lb!(net, "EX_lac__D_e", -1.0 * 60 * gf)
        lb!(net, "EX_malt_e", -0.1 * 60 * gf)
        lb!(net, "EX_gal_e", -0.2 * 60 * gf)
        lb!(net, "EX_glyc_e", -0.6 * 60 * gf)
        lb!(net, "EX_ac_e", -1.5 * 60 * gf)
        ub!(net, "EX_ac_e", 2.0 * 60 * gf)
        
        ax1 = Axis(f[i,1]; 
            title = netid,
            xlabel = "var index (sorted)", 
            ylabel = "lb shadow prices"
        )

        obj_mv = Float64[]
        obj_errv = Float64[]
        for exch in exchs
            opm = FBAOpModel(net, Clp.Optimizer)
            test_points = lb(net, exch) .* [1.0, 0.95, 0.9]
            @time obj_m, obj_err, vars_ms, vars_errs = bound_dual_prices(
                opm, exch, test_points, :lb; 
                dovars = true
            )

            sperm = sortperm(vars_ms)
            xs = eachindex(sperm)
            ys = vars_ms[sperm]
            errs = vars_errs[sperm]
            lines!(ax1, xs, ys; label = exch, linewidth = 4)
            errorbars!(ax1, xs, ys, errs)

            push!(obj_mv, obj_m)
            push!(obj_errv, obj_err)
        end
        axislegend(ax1; position = :rb)

        sperm = sortperm(obj_mv)
        ax2 = Axis(f[i,2]; 
            # aspect = AxisAspect(1),
            xticks = (eachindex(exchs), exchs[sperm]),
            xlabel = "exch", 
            ylabel = "z price", 
            xticklabelrotation = pi/4, 
            limits = (nothing, nothing, -0.23, 0.03)
        )
        barplot!(ax2, eachindex(exchs), obj_mv[sperm];
            title = "Stacked bars"
        )
    end
    f
end

## ------------------------------------------------------------