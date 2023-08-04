## ------------------------------------------------------------
@time begin
    using Plots
    using Random
    using Base.Threads
    using MetXBase
    using Statistics
    using Combinatorics
end

## ------------------------------------------------------------
struct Histogram{T}
    bin_rule::Function
    count_dict::Dict{T, Int}
    function Histogram{T}(bin_rule::Function) where T
        new{T}(bin_rule, Dict{T, Int}())
    end
end

# ------------------------------------------------------------
function find_bin(h::Histogram{T}, v) where T
    x::T = h.bin_rule(h, v) # apply bin_rule
    isa(x, T) || error("Invalid bin_rule return type, expected: ", T, ", got: ", typeof(x))
    return x
end

import Base.keys
keys(h::Histogram) = keys(h.count_dict)
import Base.values
values(h::Histogram) = values(h.count_dict)

bins(h::Histogram) = keys(h)
counts(h::Histogram) = values(h)

import Base.count!
function count!(h::Histogram{T}, v) where T
    x = find_bin(h, v)
    get!(h.count_dict, x, 0)
    h.count_dict[x] += 1
    return h
end

# Merge histograms
function count!(h0::Histogram{T}, h1::Histogram{T}, hs::Histogram{T}...) where T
    for (x, c) in h1.count_dict
        get!(h0.count_dict, x, 0)
        h0.count_dict[x] += c
    end
    for hi in hs
        for (x, c) in hi.count_dict
            get!(h0.count_dict, x, 0)
            h0.count_dict[x] += c
        end
    end
    return h0
end


## ------------------------------------------------------------
let
    kosets = []
    for it in 1:100
        push!(kosets, rand(1:100, rand(4:10)))
    end
    kosets
end

## ------------------------------------------------------------
let
    _bins = range(-1e10, 1e10; step = 0.1)
    h0 = Histogram{Float64}() do _h, v
        ci = MetXBase._find_nearest(v, _bins)
        return _bins[ci]
    end

    h_pool = [deepcopy(h0) for th in 1:4]

    T = 1 # K
    @threads for _ in 1:4
        h = h_pool[threadid()]
        for it in 1:10000000
            p = rand()
            count!(h, - log(p) * T)
        end
    end
    count!(h0, h_pool...) # reduce

    @show sum(values(h0))
    xs = sort(collect(keys(h0)))
    ws = [h0.count_dict[x] for x in xs]
    plot(xs, ws; label = "", 
        # xlim = [0,10]
    )

end

## ------------------------------------------------------------