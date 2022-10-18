## ------------------------------------------------------------------
function fixxed_mat(net::MetNets.MetNet, target_id, biom_id, bound_bins;
        rtol = 1e-2,
        fba_fun::Function = MetNets.fba, 
        solver = MetLP.Ipopt.Optimizer,  
        tocollect = eachindex(net.rxns), 
        lp_model = nothing
    )

    # idxs
    tocollect = MetNets.rxnindex.([net], tocollect)
    target_id = MetNets.rxnindex(net, target_id)
    biom_id = MetNets.rxnindex(net, biom_id)

    # flux mat (delta bin, rxns)
    Δmat = zeros(length(bound_bins), length(tocollect))

    # lp_model
    lp_model = isnothing(lp_model) ? MetLP.build_lp_model(net, solver) : lp_model

    # initial bounds
    biom_lb0, biom_ub0 = MetNets.bounds(lp_model, biom_id)
    target_lb0, target_ub0 = MetNets.bounds(lp_model, target_id)

    # fixxed loop
    lb0, ub0 = minimum(bound_bins), maximum(bound_bins)
    abstol = abs((ub0 - lb0) * rtol)

    for (i, v) in enumerate(bound_bins)
        l = max(v - (abstol/2), lb0)
        u = min(v + (abstol/2), ub0)

        MetNets.bounds!(lp_model, target_id, l, u)
        sol = fba_fun(net, biom_id; lp_model)

        fluxs = MetNets.av(sol, tocollect)
        Δmat[i, :] .= fluxs
        
    end

    # restore
    MetNets.bounds!(lp_model, biom_id, biom_lb0, biom_ub0)
    MetNets.bounds!(lp_model, target_id, target_lb0, target_ub0)

    return Δmat, lp_model
end