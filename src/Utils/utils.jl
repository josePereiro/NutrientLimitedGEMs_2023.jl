## ------------------------------------------------------------------
# Assume, All exchanges are define A -> ∅
function find_uptake_exchs(model)
    exchs = find_exchange_reaction_ids(model)
    return filter(exchs) do exch
        lb, ub = bounds(model, exch)
        return lb < 0
    end
end

## ------------------------------------------------------------------
function deep_hash(obj)
    h = hash(0.0)
    for f in fieldnames(typeof(obj))
        fv = getfield(obj, f)
        h = hash(fv, h)
    end
    return h
end

## ------------------------------------------------------------------
function format_ticks(f::Function, vals, n)
    idxs = unique!(floor.(Int, range(firstindex(vals), lastindex(vals); length = n)))
    vals = collect(vals)[idxs]
    ticks = map(f, vals)
    return (vals, ticks)
end

## ------------------------------------------------------------------
# assume lb > 0
function range_samples(net::MetNets.MetNet, rxnid, nbins::Int; 
        bound_frac = 0.9
    )
    # create range
    _, target_b0 = MetNets.bounds(net, rxnid)
    return range(target_b0 * bound_frac, target_b0; length = nbins)
end


# ------------------------------------------------------------------
function _do_th_sampling(f::Function, samples, rtol)

    lb0, ub0 = minimum(samples), maximum(samples)
    abstol = abs((ub0 - lb0) * rtol)

    @threads for (i, v) in collect(enumerate(samples))
        l = max(v - (abstol/2), lb0)
        u = min(v + (abstol/2), ub0)

        f(i, v, l, u)
    end
end

# ------------------------------------------------------------------
function _find_between(vec, i0, i1)
    @assert i0 < i1
    v0, v1 = extrema(vec)
    if i0 == 0
        v0 = zero(v0)
    elseif i0 < 0
        v0 = v0 > 0 ? zero(v0) : (v0 + abs(v0) * (1.0 + i0))
    else
        v0 = v1 < 0 ? zero(v0) : v1 - v1 * (one(i0) - i0)
    end
    
    if i1 == 0
        v1 = zero(v1)
    elseif i1 < 0
        v1 = v0 > 0 ? zero(v1) : (v0 + abs(v0) * (one(i1) + i1))
    else
        v1 = v1 < 0 ? zero(v1) : v1 - v1 * (one(i1) - i1)
    end

    @assert v0 <= v1
    return findall(vec) do v
        v >= v0 && v <= v1
    end
end