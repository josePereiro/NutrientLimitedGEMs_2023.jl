@time begin
    using Plots
    using COBREXA
    using NutrientLimitedGEMs
    const NL = NutrientLimitedGEMs
    using GLPK
    using ProjAssistant
end

## ------------------------------------------------------------------
# tools
function _normalize_mat!(mat, norm = (ubs .- lbs))
    # normalize
    for (i, row) in enumerate(eachrow(mat))
        v0 = minimum(row)
        row = abs.((row .- v0) ./ norm[i])
        mat[i, :] .= row
    end
end

# assume normalized
function _filter_noise!(mat, th)
    for (ri, row) in enumerate(eachrow(mat))
        ref = mat[ri, 1]
        for (ci, val) in enumerate(row)
            if isnan(ref) 
                ref = mat[ri, ci]
            elseif abs((ref - val) / (ref + val)) > th
                mat[ri, ci] = NaN
            else
                ref = mat[ri, ci]
            end
        end
    end
end

# assume normalized
function _round_zeros!(mat, th)
    for (ri, row) in enumerate(eachrow(mat))
        for (ci, val) in enumerate(row)
            if abs(mat[ri, ci]) <= th
                mat[ri, ci] = zero(th)
            end
        end
    end
end

## ------------------------------------------------------------------
# Heatmap
let
    return 
    
    # load globals
    model_params, nsample = ldat(NL, "effect_mat-globals.jls")
    
    # load model
    @time model = load_HartGEM_model(; model_params...)
    lbs, ubs = bounds(model)

    # compute dep mat
    for exchid in find_uptake_exchs(model)
        println("\n", "-"^60)
        
        @time vals, mat = load_compute_effect_mat_cache(model, exchid)

        # format
        norm = (ubs .- lbs)
        _normalize_mat!(mat, norm)
        _filter_noise!(mat, 0.4)
        _round_zeros!(mat, 1e-5)

        # pretty
        names = lrawdat(NL, "met_readable_ids.jls")
        format_id = (met) -> get(names, met, met)

        # plot
        p = heatmap(eachindex(vals), eachindex(model.rxns), mat; 
            xticks = format_ticks(i -> round(vals[i]; sigdigits = 2), eachindex(vals), 4),
            xrotation = 35,
            title = string(exchid), 
            xlabel = string("flux value: ", reaction_eq(model, exchid; format_id)), 
            ylabel = "all reactons"
        )

        sfig(NL, p, "preview", ".png")
        sfig(NL, p, 
            @fileid, "HartGEM-flux-dep-matrix", model_params, (;exchid), ".png";
            verbose = true
        )
    end
end

## ------------------------------------------------------------------
# _exchs_fixxed_trend
function _exchs_fixxed_trend(model, biomider, model_params; 
        eps = 1e-5, 
        plot_biom = true,
        plot_exchs = true
    )

    biomidx = rxnindex(model, biomider)
    lbs, ubs = bounds(model)
    M, N = size(model.S)

    # compute dep mat
    exch_ids = find_uptake_exchs(model)
    for exchid in exch_ids

        println("\n", "-"^60)
        println(exchid)
        
        @time _, mat = load_compute_effect_mat_cache(model, exchid)
        exch_vals = mat[rxnindex(model, exchid), :]

        # format
        # norm = fill((ubs[biomidx] .- lbs[biomidx]), N)
        # _normalize_mat!(mat, norm)
        # _filter_noise!(mat, 0.4)
        # _round_zeros!(mat, 1e-5)
        
        # exchanges

        # pretty
        names = lrawdat(NL, "met_readable_ids.jls")
        format_id = (met) -> get(names, met, met)

        # PLOT
        p = plot()
        # exchanges
        if plot_exchs
            for id in exch_ids
                idx = rxnindex(model, id)
                # exch_vec = log10.(mat[idx, :] .+ eps)
                exch_vec = mat[idx, :]
                exch_vec .= exch_vec ./ maximum(abs, exch_vec)
                plot!(p, exch_vals, exch_vec; label = id, alpha = 0.8, lw = 3)
            end
        end

        # biom
        if plot_biom
            # biom_vec =  log10.(mat[biomidx, :] .+ eps)
            biom_vec =  mat[biomidx, :]
            biom_vec .= biom_vec ./ maximum(abs, biom_vec)
            plot!(p, exch_vals, biom_vec; lw = 4, c = :white, label = "")
            plot!(p, exch_vals, biom_vec; lw = 2, c = :red, label = "biom")
        end

        plot!(p;
            # xticks = format_ticks(i -> round(vals[i]; sigdigits = 2), eachindex(vals), 4),
            xlabel = string("flux value (", reaction_eq(model, exchid; format_id), ")"), 
            # ylabel = "normalized value (all reactions)",
            ylabel = "reaction val",
            title = string(exchid), 
        )

        # sfig(NL, p, "preview", ".png")
        sfig(NL, p, 
            @fileid, "flux-dep-trend", model_params, (;exchid), ".png";
            verbose = true
        )
    end
end

## ------------------------------------------------------------------
# biom_trend
function _exchs_biom_trend(model, biomider, model_params)
    
    biomidx = rxnindex(model, biomider)

    # compute dep mat
    for exchid in find_uptake_exchs(model)

        println("\n", "-"^60)
        println(exchid)
        
        @time vals, mat = load_compute_effect_mat_cache(model, exchid)
        exchidx = rxnindex(model, exchid)
        exch_vals = mat[exchidx, :]

        # format
        mat = round.(mat; digits = 3)

        # pretty
        names = lrawdat(NL, "met_readable_ids.jls")
        format_id = (met) -> get(names, met, met)

        # PLOT
        p = plot()

        # biom
        biom_mat = mat[biomidx, :]
        plot!(p, exch_vals, biom_mat; lw = 4, c = :white, label = "")
        plot!(p, exch_vals, biom_mat; lw = 2, c = :red, label = "biom")
        scatter!(p, exch_vals, biom_mat; ms = 4, c = :red, label = "")
        
        # zero marker
        plot!(p, [0, 0], [0, maximum(biom_mat)];
            label = "", c = :black, alpha = 0.5, lw = 3, ls = :dash
        )

        plot!(p;
            xlabel = string("flux value (", reaction_eq(model, exchid; format_id), ")"), 
            ylabel = "normalized value (all reactions)",
            title = string(exchid), ylim = [0, Inf], 
            legend = :bottomleft
        )

        sfig(NL, p, 
            @fileid, "biom-fix-exch-trend", model_params, (;exchid), ".png";
            verbose = true
        )
    end

end

## ------------------------------------------------------------------
# HartGEMs
let
    return

    # load globals
    globals_file = "effect_mat-globals-HartGEMs.jls"
    model_params, nsamples = ldat(NL, globals_file)

    # load model
    @time model = load_HartGEM_model(; model_params...)
    biomider = HUMAN_BIOMASS_IDER

    _exchs_biom_trend(model, biomider, model_params)

end

## ------------------------------------------------------------------
# EColi
let
    return

    # load globals
    globals_file = "effect_mat-globals-EColi.jls"
    model_params, nsamples = ldat(NL, globals_file)

    # load model
    biomider = "R_BIOMASS_Ecoli_core_w_GAM"
    @time model = load_ecoli()

    _exchs_biom_trend(model, biomider, model_params)
end


## ------------------------------------------------------------------
# ToyModel
let
    # load globals
    globals_file = "effect_mat-globals-ToyModel.jls"
    model_params, nsamples = ldat(NL, globals_file)

    # load model
    biomider = "z"
    @time model = NL.load_toy_model3()

    # _exchs_biom_trend(model, biomider, model_params)
    _exchs_fixxed_trend(model, biomider, model_params)
end
