## ------------------------------------------------------------------
@time begin
    using NutrientLimitedGEMs_2023
    const NL = NutrientLimitedGEMs_2023

    using ProjFlows
    using Gurobi
    using MetXBase
    
    using Plots
    
end

# ------------------------------------------------------------------
include("0_params.jl")
include("1_utils.jl")

## ------------------------------------------------------------------
let

    D = 0.5 # 1/ h

    X_ts = Float64[]

    c_glc0 = 12.0
    c_nh40 = 8.0
    c_glc_ts = [
        ones(20); 
        1.0 .- sin.((0:52) .* 6e-1) .* 0.08
        ones(20); 
        ones(50); 
        ones(20);
    ] .* c_glc0
    c_nh4_ts = [
        ones(20); 
        ones(50);
        ones(20);
        1.0 .- sin.((0:52) .* 6e-1) .* 0.12;
        ones(20);
    ] .* c_nh40

    ps = Plots.Plot[]

    for (i, (c_glc, c_nh4)) in enumerate(zip(c_glc_ts, c_nh4_ts))

        net0 = pull_net("ecoli_core")
        lep0 = lepmodel(net0)
        biom_id = extras(lep0, "BIOM")
        glc_id = extras(lep0, "EX_GLC")
        nh4_id = extras(lep0, "EX_NH4")
        bounds!(lep0, biom_id, D * 0.99, D * 1.01)
        
        opm = FBAOpModel(lep0, LP_SOLVER)
        
        set_linear_obj!(opm, glc_id, -1.0)
        optimize!(opm)
        max_q_glc = solution(opm, glc_id)
        
        set_linear_obj!(opm, nh4_id, -1.0)
        optimize!(opm)
        max_q_nh4 = solution(opm, nh4_id)

        X1 = abs(c_glc * D / max_q_glc)
        X2 = abs(c_nh4 * D / max_q_nh4)
        
        X = min(X1, X2)
        # @show X
        push!(X_ts, X)

        p1 = plot(; xlabel = "time", ylabel = "biomass")
        plot!(p1, X_ts;
            label = "",
            ylim = (0.0, 1.0), 
            c = :black,
            lw = 3
        )
        
        p2 = plot(; xlabel = "time", ylabel = "fresh medium conc")
        plot!(p2, c_glc_ts[1:i]; 
            label = "c_glc", 
            ylim = (0.0, 20.0),
            c = :red, 
            lw = 3
        )
        plot!(p2, c_nh4_ts[1:i]; 
            label = "c_nh4", 
            ylim = (0.0, 20.0),
            c = :blue, 
            lw = 3
        )
        
        p = plot(p2, p1; 
            size = (800, 400), 
            margin = 5.0Plots.mm
        )
        
        push!(ps, p)
    end

    sgif(PROJ, ps, "nut_lim.gif") |> println
    sfig(PROJ, ps[end], "nut_lim.png") |> println

end