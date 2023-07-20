## ------------------------------------------------------------------
@time begin
    using ProjFlows
    using MetXBase
    using Plots
    using Random
    using Distributions
    using Base.Threads
    using ForwardDiff
    using Optim
end

# ------------------------------------------------------------------
include("0_params.jl")
include("1_utils.jl")

## ------------------------------------------------------------------
# Multivariate
let
    # Define the function to minimize
    # Random.seed!(123)
    M = rand(5:20)
    D = rand(5:10)
    global _mix = _rand_mixture(M, D; 
        μscale = 1.0, σ0 = 0.3,
    )

    # solver = GradientDescent()
    # solver = BFGS()
    solver = LBFGS()
    # solver = ConjugateGradient()
    # solver = Newton()
    # solver = NelderMead()
    @time global _max_mix = _max_mixture_H(_mix, _entropy_lb1, solver)
    _suggest_range(_mix, 2)
    
    @show _entropy_lb1(_mix)
    @show _entropy_lb1(_max_mix)
    @show _entropy_mcf(_mix, Int(1e5))
    @show _entropy_mcf(_max_mix, Int(1e5))
    @show _entropy_ub1(_mix)
    @show _entropy_ub1(_max_mix)
    

end

## ------------------------------------------------------------------
## ------------------------------------------------------------------
# TODO: redo the validation script using Optim
# TODO: make the MvNormal version
# TODO: Think about the truncated version
let

    # Define the function to minimize
    Random.seed!(123)
    M = 15
    qs_t = [Normal(5.0 * randn(), 0.1 + rand()) for m in 1:M]
    uni_mix = MixtureModel(qs_t)
    
    # solver = GradientDescent()
    # solver = BFGS()
    solver = LBFGS()
    # solver = ConjugateGradient()
    # solver = Newton()
    # solver = NelderMead()
    @time maxHl_mix = _max_mixture_H(uni_mix, _entropy_lb1, solver)
    @time maxHu_mix = _max_mixture_H(uni_mix, _entropy_ub1, solver)

    
    for (id, md) in [
            ("uni_mix", uni_mix),
            ("maxHl_mix", maxHl_mix),
            ("maxHl_mix", maxHu_mix),
        ]
        println("-"^40)
        @info(id,
            Hl = _entropy_lb1(md),
            Hmu = _entropy_mcu(md, Int(1e5)),
            Hmf = _entropy_mcf(md, Int(1e5)),
            Hq = _entropy_quad(md),
            Hu = _entropy_ub1(md),
        )
    end

    
    # @show maxH_mix.ws

    p = plot()
    _plot!(p, uni_mix; lw = 2, c = :black, label = "uniform")
    _plot!(p, maxHl_mix; lw = 2, c = :blue, label = "max Hl")
    _plot!(p, maxHu_mix; lw = 2, c = :red, label = "max Hu")
    p

end


## ------------------------------------------------------------------
let
    # Define the function to minimize
    # Random.seed!(123)
    M = 3
    qs_t = [Normal(1.0 * randn(), rand()) for m in 1:M]
    uni_mix = MixtureModel(qs_t)
    
    n = 10^6
    Hmu = _entropy_mcu(uni_mix, n)
    Hmf = _entropy_mcf(uni_mix, n)

    Hl = _entropy_lb1(uni_mix)
    Hq = _entropy_quad(uni_mix)
    Hu = _entropy_ub1(uni_mix)
    
    @show Hl Hq Hmu Hmf Hu
    
    p =plot()
    _plot!(p, uni_mix)


    # @time result = _monte_carlo_integration_lazy(f, a, b, n)
    # println("The result of the integration is: ", result)
end

## ------------------------------------------------------------------
## ------------------------------------------------------------------
let
    Random.seed!(123)

    M = 3
    
    qs_t = [Normal(5.0 * randn(), 0.3 + rand()) for m in 1:M]
    ws = rand(M)
    ws ./= sum(ws)
    @show ws
    
    f = (_ws) -> let 
        nm = MixtureModel(_ws, qs_t)
        return _entropy_lb1(nm)
        # sum(_ws .* ws)
    end

    @time ws1 = ForwardDiff.gradient(f, ws)
    @time ws1 = ForwardDiff.gradient(f, ws)
    @show ws1

    nothing
end


# @show Hl = _entropy_lb1(nm)
# @show Hq = _entropy_quad(nm, Int(1e4))
# @show Hu = _entropy_ub1(nm)