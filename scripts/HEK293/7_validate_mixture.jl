## ------------------------------------------------------------------
@time begin
    using ProjFlows
    using MetXBase
    using Plots
    using Random
    using Distributions
    using Base.Threads
end

# ------------------------------------------------------------------
include("0_params.jl")
include("1_utils.jl")

## ------------------------------------------------------------------
# mc vd D
let
    return 
    Random.seed!(321)
    M = 10
    ns = round.(Int, 10.0.^[1:0.5:6;])
    p = plot()
    for D in [1, 3, 5, 10]
        Hmf = zeros(length(ns))
        nm = _rand_mixture(M, D; 
            μscale = 1.0, σ0 = 0.3,
        )
        @threads for (i, n) in collect(enumerate(ns))
            # D = rand(5:10)
            Hmf[i] = _entropy_mcf(nm, n)
        end

        plot!(p, log10.(ns), Hmf .- Hmf[1]; 
            label = string("D: ", D), lw = 2
        )
    end
    p
end

## ------------------------------------------------------------------
# Correlation Hq vs Hl and Hu
let
    Random.seed!(123)
    NT = Int(1e3)
    for D in [50, 100, 500]
        Hls = zeros(NT)
        Hmf = zeros(NT)
        Hus = zeros(NT)
        @threads for t in 1:NT
            M = rand(5:20)
            nm = _rand_mixture(M, D;
                μscale = 1.0, σ0 = 0.3,
            )
            Hls[t] = _entropy_lb1(nm)
            # TODO: explore at which n _entropy_mcf saturate with respect to 
            Hmf[t] = _entropy_mcf(nm, Int(1e4))
            Hus[t] = _entropy_ub1(nm)
        end

        p = plot(;xlabel = "stand Hmc", ylabel = "stand Hl | Hu")
        _Hqs = _rescale(Hmf)
        
        _Hus = _rescale(Hus)
        scatter!(p, _Hqs, _Hus; 
            label = string("Hu (ρ=", round(cor(_Hqs, _Hus); sigdigits = 2), ")"), 
            c = :red,
            msc=:auto,
            alpha = 0.5
        )

        _Hls = _rescale(Hls)
        scatter!(p, _rescale(Hmf), _rescale(Hls); 
            label = string("Hl (ρ=", round(cor(_Hqs, _Hls); sigdigits = 2), ")"), 
            c = :blue,
            msc=:auto,
            alpha = 0.5
        )

        plot!(p, sort(_Hqs), sort(_Hqs); 
            label = "Hmc = Hmc", lw = 3, 
            ls = :dash, c = :black, 
            alpha = 0.7
        )
        
        sfig(PROJ, p, 
            ["Mixture"], "Hl_Hu_Hmc_correlation", (;D), ".png"
        ) |> println
    end
end

exit()

## ------------------------------------------------------------------
let 
    M = 5
    qs = [Normal(randn(), rand()) for m in 1:M]
    
    NT = Int(1e5)

    # MC
    Uniform_nm = NormalMixture(qs)
    H = _entropy_lb1(Uniform_nm)
    maxH_ws = Uniform_nm.ws
    lk = ReentrantLock()
    @time @threads for t in 1:NT
        ws_t = rand(M)
        ws_t ./= sum(ws_t)
        nm = NormalMixture(ws_t, qs)
        Ht = _entropy_lb1(nm)
        if Ht > H
            lock(lk) do
                H = Ht
                maxH_ws = ws_t
                @show H
            end
        end
    end

    maxH_nm = NormalMixture(maxH_ws, qs)
    
    @show _entropy_lb1(maxH_nm)
    @show _entropy_quad(maxH_nm)
    @show _entropy_lb1(Uniform_nm)
    @show _entropy_quad(Uniform_nm)

    ps = Plots.Plot[]
    p = plot()
    @time _plot!(p, Uniform_nm; 
        label = "uniform", c = :blue, lw = 3, 
    )
    # _plot!.([p], Uniform_nm.qs; 
    #     label = "", c = :blue, lw = 3, 
    #     alpha = 0.5
    # )
    push!(ps, p)

    p = plot()
    @time _plot!(p, maxH_nm; 
        label = "max Hl", c = :red, lw = 3, 
    )
    # _plot!.([p], maxH_nm.qs; 
    #     label = "", c = :red, lw = 3, 
    #     alpha = 0.5
    # )
    push!(ps, p)
    
    
    plot(ps...)
end


## ------------------------------------------------------------------
## ------------------------------------------------------------------
# Gradiend descent
let
    Random.seed!(190)

    ## ---------------------------------------------
    # Mixure
    M = 50
    qs = [Normal(randn() * 2.0, 0.2 + rand()) for m in 1:M]

    ## ---------------------------------------------
    # GD
    target = ones(M) .* 0.0
    x0 = rand(M)
    x0 ./= sum(x0)
    x1 = rand(M)
    x1 ./= sum(x1)
    gdth = 1e-5
    maxiter = Int(1e5)
    verbose = false
    maxΔx = ones(M) .* 1e-6
    minΔx = ones(M) .* 1e-7

    function up_fun(gdmodel) 
        # Up weight
        _ws = gd_value(gdmodel)
        _ws .= abs.(_ws)
        _ws ./= sum(_ws)

        _maxHl_nm = NormalMixture(_ws, qs)
        # Hl = _entropy_lb1(_maxHl_nm)
        Hl = _entropy_quad(_maxHl_nm)
        
        return fill(Hl, M)
    end

    gdmodel = grad_desc_vec(up_fun; 
        target, x0, x1, minΔx, maxΔx, gdth, maxiter, verbose
    )

    # Hl = up_fun(gdmodel)

    maxHl_nm = NormalMixture(gd_value(gdmodel), qs)
    Uniform_nm = NormalMixture(qs)

    @show _entropy_lb1(maxHl_nm)
    @show _entropy_lb1(Uniform_nm)
    @show _entropy_quad(maxHl_nm)
    @show _entropy_quad(Uniform_nm)

    ps = Plots.Plot[]
    p = plot()
    _plot!(p, maxHl_nm; label = "max Hl", lw = 2, c = :red)
    _plot!(p, Uniform_nm; label = "Uniform_nm Hl", lw = 2, c = :blue)
    push!(ps, p)
    
    p = plot()
    bar!(p, sort(maxHl_nm.ws); label = "max Hl ws", lw = 2, c = :red)
    plot!(p, Uniform_nm.ws; label = "Uniform_nm ws", lw = 2, c = :blue)
    push!(ps, p)
    
    plot(ps...)

end

## ------------------------------------------------------------------
## ------------------------------------------------------------------
let 
    M = 3
    qs = [Normal(randn(), rand()) for m in 1:M]
    global nm = NormalMixture(qs)

    @show _entropy_lb1(nm)
    p = plot()
    foreach(nm.qs) do q
        _plot!(p, q; label = "")
    end
    _plot!(p, nm; 
        label = "mixture", lw = 3, alpha = 0.4
    )
    p
    
end

"""
Here is a long poem with regular metric based on the title and abstract of the research paper:

**Inference of metabolic fluxes in nutrient-limited continuous cultures**

A challenge in the study of life
Is to infer the metabolic strife
That cells endure in different states
With limited experimental traits

How can we make predictions sound
From measurements that are not abound?
How can we use the data we have
To understand the metabolic salve?

One way to do this is to apply
The principle of Maximum Entropy
Which says that from a set of constraints
We choose the most likely and consistent states

This method has been used before
To model metabolism at its core
But can we make it even better
By using less data and more clever?

We propose a thesis that is bold
That chemostat dynamics can unfold
The secrets of the nutrient-limited cells
By driving them to states that are easy to tell

We say that in these cultures we can see
That growth rate and uptake rates are the key
To infer the fluxes that we need
Without measuring them indeed

We test our model with simulations first
And then with real data we quench our thirst
We compare our results with other ways
That use convex optimization in their phase

We find that our approach is quite robust
And can infer the fluxes with more trust
We also show that it is more efficient
And can use less data with more intent

We hope that our work can inspire more
To use Maximum Entropy at their core
And to explore the chemostat dynamics more
To reveal the metabolic secrets galore.
"""