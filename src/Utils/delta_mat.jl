# ------------------------------------------------------------------
function delta_mat(
        net::MetNets.MetNet, 
        target_id, biom_id, bound_bins;
        fba_fun::Function = MetLP.r2_fba, 
        solver = MetLP.Ipopt.Optimizer, 
        tocollect = eachindex(net.rxns), 
        lp_model = nothing
    )

    # TODO: assert all(lb .== 0.0)

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

    # bound loop
    for (i, ub) in enumerate(bound_bins)

        MetNets.bounds!(lp_model, target_id; lb = 0.0, ub)
        sol = fba_fun(net, biom_id; lp_model)
        
        fluxs = MetNets.av(sol, tocollect)
        Δmat[i, :] .= fluxs

    end

    # restore
    MetNets.bounds!(lp_model, biom_id, biom_lb0, biom_ub0)
    MetNets.bounds!(lp_model, target_id, target_lb0, target_ub0)

    return Δmat, lp_model
end