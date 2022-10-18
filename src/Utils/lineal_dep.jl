function linear_dep(delta_mat::Matrix, ind_vals)
    N = size(delta_mat, 2)
    ms, errs = zeros(N), zeros(N)
    for (i, rxn_vals) in enumerate(eachcol(delta_mat))
        _, ms[i], errs[i] = linear_fit(ind_vals, rxn_vals)
    end
    return ms, errs
end

function linear_dep(delta_mat::Matrix, target_idx::Int)
    ind_vals = delta_mat[:, target_idx]
    return linear_dep(delta_mat, ind_vals)
end

function linear_dep(
        net::MetNets.MetNet, 
        target_id, biom_id, bound_bins; 
        mat_fun::Function = delta_mat,
        mat_kwargs...
    )

    target_id = MetNets.rxnindex(net, target_id)
    delta_mat, lp_model = mat_fun(net, target_id, biom_id, bound_bins; 
        mat_kwargs...
    )
    ms, errs = linear_dep(delta_mat, bound_bins)
    return ms, errs, delta_mat, lp_model
end

function linear_dep(
        net::MetNets.MetNet, 
        target_ids::Vector, biom_id;
        nbins::Int = 5, 
        bound_frac::Float64 = 0.9,
        lp_model = nothing,
        mat_kwargs...
    )

    ms_mat, errs_mat = nothing, nothing
    for (i, tid) in enumerate(target_ids)
        
        # create range
        _, target_b0 = MetNets.bounds(net, tid)
        bound_bins = range(target_b0 * bound_frac, target_b0; length = nbins)

        ms, errs, _, lp_model = linear_dep(net, tid, biom_id, bound_bins; 
            lp_model, mat_kwargs...
        )

        # init holder
        if isnothing(ms_mat)
            ms_mat = zeros(length(ms), length(target_ids))
            errs_mat = zeros(length(ms), length(target_ids))
        end
        ms_mat[:, i] .= ms
        errs_mat[:, i] .= errs
    end

    return ms_mat, errs_mat, lp_model

end