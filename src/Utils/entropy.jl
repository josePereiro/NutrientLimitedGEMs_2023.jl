function uniform_entropy(
        net::MetNets.MetNet;
        alpha = Inf, 
        clear_cache = false, 
        do_fva = true,
        tocollect = [],
        verbose = false,
        onerr::Function = rethrow,
        fva_eps = 1e-2,
        cache_seed = :uniform_entropy,
        converge_ep_kwargs...
    )

    net_hash = deep_hash(net)

    try

        # fva_preprocess
        # cid = (cache_seed, net_hash, :fva_preprocess, fva_eps)
        # clear_cache && delcache(NutrientLimitedGEMs, cid)
        # if do_fva
            # net = lcache(NutrientLimitedGEMs, cid) do
        net = MetLP.fva_preprocess(net; verbose, eps = fva_eps)
        #     end
        # end
        
        # EP
        # cid = (:converge_ep!, cache_seed, net_hash, converge_ep_kwargs, tocollect)
        # clear_cache && delcache(NutrientLimitedGEMs, cid)
        # S, flxs = lcache(NutrientLimitedGEMs, cid) do
            
            epmodel = MetEP.EPModel(net; alpha)
            sol = MetEP.converge_ep!(epmodel; 
                verbose, converge_ep_kwargs...
            )
            
            Q = MetEP.get_join(epmodel)
            S_ = entropy(Q)

            flxs_ = MetNets.av.([net], [sol], tocollect)

        #     return S_, flxs_
        # end
        S, flxs = S_, flxs_

        verbose && println()
        verbose && @info("Done!!!", S)
        
        return net, S, flxs

    catch err
        (err isa InterruptException) && rethrow(err)
        onerr(err)
        @warn err
        return net, NaN, NaN
    end

end


function delta_entropy(
        net::MetNets.MetNet, target_id, bound_bins;
        verbose = false,
        cache_seed = :delta_entropy,
        uniform_entropy_kwargs...
    )

    net_hash = deep_hash(net)
    Svec = nothing
    flxs_mat = nothing

    for (i, ubi) in enumerate(bound_bins)

        verbose && println("\n", "="^60)
        verbose && @info("Doing", ub = ubi)
        verbose && println()

        # local net
        net_ = MetNets.MetNet(net; 
            lb = deepcopy(net.lb),
            ub = deepcopy(net.ub),
        )
        MetNets.bounds!(net_, target_id; lb = 0.0, ub = ubi)

        S, flxs = uniform_entropy(net_; 
            verbose,
            cache_seed = (cache_seed, net_hash),
            uniform_entropy_kwargs...
        )

        if isnothing(flxs_mat)
            Svec = zeros(length(bound_bins))
            flxs_mat = zeros(length(bound_bins), length(flxs))
        end
        
        Svec[i] = S
        flxs_mat[i, :] .= flxs

    end

    return Svec, flxs_mat
end