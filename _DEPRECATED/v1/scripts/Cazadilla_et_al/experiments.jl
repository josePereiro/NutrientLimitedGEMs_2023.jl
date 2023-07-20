# TODO: Add to MetXCultureHub

# ------------------------------------------------------------------
@time begin
    using MetX
    using NutrientLimitedGEMs
    using CSV
    using DataFrames
    using ProjFlows
    using ProjFlows: rand_color, rand_line_style, rand_marker_shape
    using Plots
    using Statistics
end

# ------------------------------------------------------------------
PROJ = Project0(NutrientLimitedGEMs)

# ------------------------------------------------------------------
let
    fn = "/Users/Pereiro/.julia/dev/NutrientLimitedGEMs/data/processed/Calzadilla_et_al/Extracted.csv"
    # Filename,D(1/h),Fecha,hora,Tiempo cultivo (d),Xv,ALA:[area],ALA:[C(mM)],ALA:[qS(Area/cel*h)],ARG:[area],ARG:[C(mM)],ARG:[qS(Area/cel*h)],ASN:[area],ASN:[C(mM)],ASN:[qS(Area/cel*h)],ASP:[area],ASP:[C(mM)],ASP:[qS(Area/cel*h)],GLN:[area],GLN:[concentracion],GLN:[qS(nmol/e6cel*h)],GLN:[qS(area)],GLN:[DCo/Xv],GLU:[area],GLU:[C(mM)],GLU:[qS(Area/cel*h)],GLC:[area],GLC:[concentracion],GLC:[qS(nmol/e6cel*h)],GLC:[qS(area)],GLC:[DCo/Xv],GLY:[area],GLY:[C(mM)],GLY:[qS(Area/cel*h)],ILE:[area],ILE:[C(mM)],ILE:[qS(Area/cel*h)],LAC:[area],LAC:[concentracion],LAC:[qS(nmol/e6cel*h)],LAC:[qS(area)],LEU:[area],LEU:[C(mM)],LEU:[qS(Area/cel*h)],LYS:[area],LYS:[C(mM)],LYS:[qS(Area/cel*h)],MET:[area],MET:[C(mM)],MET:[qS(Area/cel*h)],PHE:[area],PHE:[C(mM)],PHE:[qS(Area/cel*h)],PRO:[area],PRO:[C(mM)],PRO:[qS(Area/cel*h)],SER:[area],SER:[C(mM)],SER:[qS(Area/cel*h)],THR:[area],THR:[C(mM)],THR:[qS(Area/cel*h)],TRP:[area],TRP:[C(mM)],TRP:[qS(Area/cel*h)],TYR:[area],TYR:[C(mM)],TYR:[qS(Area/cel*h)],VAL:[area],VAL:[C(x)],VAL:[qS]
    global df = CSV.read(fn, DataFrame)
    
    # Utils
    global metids = ["ALA", "ARG", "ASN", "ASP", "GLN", "GLU", "GLC", "GLY", "ILE", "LAC", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
    global met_colors_vec = Plots.distinguishable_colors(length(metids))
    global met_colors_dict = Dict((metids .=> met_colors_vec)...)
    
    df[!, "StSt"] = zeros(Int, size(df, 1))
    df[7:12, "StSt"] .= 1
    df[26:29, "StSt"] .= 2
    df[34:38, "StSt"] .= 3
    df[47:54, "StSt"] .= 4
    df[58:61, "StSt"] .= 5
    
    # Convert
    df[!, "D[1/d]"] .= df[!, "D[1/h]"] .* 24
    
    # Cell biomass
    # Xv [cell/ml] * 1e3 = Xv [cell/ L]
    # Xv [cell/L] * CDW [pg/ cell] = Xv [pg/L]
    # Xv [pg/L] * 1e-12 = Xv [gDW/L]
    # CDW = 514 # pg/cell from 10.1371/journal.pone.0043394
    CDW = 230 # [pg/cell] To improve max biom
    df[!, "Xv[gDW/L]"] = df[:, "Xv[cell/ml]"] .* 1e3 .* CDW .* 1e-12

    # Fluxes
    # q = (c - s) D / Xv [mmol/L × h^{-1} × L/gDW]
    for id in metids
        c = df[1, "$id:C[mM]"]
        D = df[!, "D[1/h]"]
        Xv = df[!, "Xv[gDW/L]"]
        ss = df[!, "$id:C[mM]"]
        df[!, "$id:q[mmol/gDW h]"] = -(c .- ss) .* D ./ Xv
    end

    # StSt bundle
    global StSt_dat = Dict()
    
    # Utils
    StSt_dat["metids"] = metids
    StSt_dat["StSts"] = 1:5

    for st in StSt_dat["StSts"]
        
        dat = Dict()
        
        # Utils
        dat["metids"] = metids
        idxs = dat["idxs"] = findall(df[!, "StSt"] .== st)
        
        # D [1/h]
        dat["D"] = Dict(
            "unit" => "h^{-1}",
            "val" => mean(df[idxs, "D[1/h]"]),
            "err" => std(df[idxs, "D[1/h]"])
        )

        # Xv [1/h]
        dat["Xv"] = Dict(
            "unit" => "gDW L^{-1}",
            "val" => mean(df[idxs, "Xv[gDW/L]"]),
            "err" => std(df[idxs, "Xv[gDW/L]"])
        )

        for id in metids

            # c [mM]
            dat["c$id"] = Dict(
                "unit" => "mM",
                "val" => df[1, "$id:C[mM]"],
                "err" => 0.0
            )

            # s [mM]
            dat["s$id"] = Dict(
                "unit" => "mM",
                "val" => mean(df[idxs, "$id:C[mM]"]),
                "err" => std(df[idxs, "$id:C[mM]"])
            )

            # q [mM]
            dat["q$id"] = Dict(
                "unit" => "mM",
                "val" => mean(df[idxs, "$id:q[mmol/gDW h]"]),
                "err" => std(df[idxs, "$id:q[mmol/gDW h]"])
            )

        end

        StSt_dat[st] = dat
    end

    # Save
    sprocdat(PROJ, (;StSt_dat, df), 
        ["Calzadilla_et_al"], "Calzadilla.bundle", ".jls"; 
        verbose = true
    )
    nothing
end

# --------------------------------------------------------
# Plots

# --------------------------------------------------------
# medium
let
    cs = [df[1, "$id:C[mM]"] for id in metids]
    sidx = sortperm(cs)

    p = bar(
        eachindex(metids), cs[sidx];
        label = "", 
        xlabel = "nutrinet",
        ylabel = "medium conc [mM]",
        color = met_colors_vec[sidx], 
        xticks = (eachindex(metids), metids[sidx]), 
        xrotation = 45
    )

       # save
    sfig(PROJ, p, 
        ["Calzadilla_et_al"], "medium_concs", ".png";
    ) |> println
end

# --------------------------------------------------------
# relative concentrtion
let
    p = plot(; xlabel = "time [d]", ylabel = "log10(s/c)")
        
    for z in 1:2
        for id in metids

            xs = df[!, "Tiempo cultivo [d]"]
            ys = log10.(df[!, "$id:C[mM]"] ./ df[1, "$id:C[mM]"])

            # highlite ststs
            z == 1 && for st in 1:5
                idx = findall(df[!, "StSt"] .== st)
                isempty(idx) && continue
                x0, x1 = extrema(idx)
                vspan!(p, xs[[x0, x1]], color = :gray, fillalpha = 0.03, labels = "");
            end

            # rel concentration
            z == 2 && plot!(p, xs, ys; 
                legend=:outertopright,
                lw = 3,
                color = met_colors_dict[id],
                label = id
            )
        end
    end # z

    # save
    sfig(PROJ, p, 
        ["Calzadilla_et_al"], "Rel_conc", ".png";
    ) |> println
end

# --------------------------------------------------------
# Abs Fluxes
let
    p = plot(; xlabel = "time [d]", ylabel = "log10 abs q [mmol/ gDW h]")
    
    for z in 1:2
        for id in metids

            xs = df[!, "Tiempo cultivo [d]"]
            ys = log10.(abs.(df[!, "$id:q[mmol/gDW h]"]))

            # highlite ststs
            z == 1 && for st in 1:5
                idx = findall(df[!, "StSt"] .== st)
                isempty(idx) && continue
                x0, x1 = extrema(idx)
                vspan!(p, xs[[x0, x1]], color = :gray, fillalpha = 0.03, labels = "");
            end

            # rel concentration
            z == 2 && plot!(p, xs, ys; 
                legend=:outertopright,
                color = met_colors_dict[id],
                lw = 3,
                label = id
            )
        end
    end # z

    # save
    sfig(PROJ, p, 
        ["Calzadilla_et_al"], "Abs_fluxes", ".png";
    ) |> println
end

# --------------------------------------------------------
# StSts
let
    
    ps = Plots.Plot[]
    for id in StSt_dat["metids"]    
        vals = [StSt_dat[st]["q$id"]["val"] for st in StSt_dat["StSts"]]
        errs = [StSt_dat[st]["q$id"]["err"] for st in StSt_dat["StSts"]]
        
        p = scatter(string.("SS", StSt_dat["StSts"]), vals; 
            yerr = errs,
            title = string("q", id, "\n[mmol/ gDW h]"), 
            ylabel = "",
            label = "", 
            color = met_colors_dict[id],
            m = 9,
            thickness_scaling = 1.3,
            tickfontsize = 16,
            titlefontsize = 16,
            guidefontsize = 16,
            xrotation = 45
        )
        push!(ps, p)
    end

    # save
    sfig(PROJ, ps, 
        ["Calzadilla_et_al"], "StSts", "fluxes", ".png"
    ) |> println
end