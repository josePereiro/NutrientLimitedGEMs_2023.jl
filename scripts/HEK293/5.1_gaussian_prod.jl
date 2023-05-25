## ------------------------------------------------------------------
@time begin
    using NutrientLimitedGEMs
    const NL = NutrientLimitedGEMs

    using ProjFlows
    using MetX
    using Gurobi
    using DataFrames
    using Plots
    using Distributions
end

# -------------------------------------------------------------------
include("0_params.jl")
include("1_utils.jl")

## ------------------------------------------------------------------
@time let

    c1 = 5500
    c = 0
    biom1_th = 0.3
    traj_dir = procdir(PROJ, ["HEK293", "sims"])
    global nl_leps = LEPModel[]

    global fva_bins = Set{BitVector}()

    for fn in readdir(traj_dir; join = true)
        endswith(fn, ".jls") || continue
        global _, sim = ldat(fn)
        sim["status"] == :success || continue
        haskey(sim, SIM_ID) || continue
        
        # biom
        biom1 = sim["biom1"]
        biom1 > biom1_th || continue
        
        println("-"^60)
        @show fn
        
        # box network
        blep = sim["fva_lep"]
        
        bounded = blep.lb .< 1e-4
        bounded = bounded .& (blep.ub .< 1e-4)
        (bounded in fva_bins) && continue
        
        push!(fva_bins, bounded)
        push!(nl_leps, blep)
        
        c += 1
        c == c1 && break
    end

    global lep0 = sim["lep0"]

    nothing
end

## ------------------------------------------------------------------
# Components
@time let
    
    global Qs = MvNormal[]
    
    for lep in nl_leps
        # EP
        _, epm = withcachedat(PROJ, :get!, (:FluxEPModelT0, lep,)) do 
            epm_ = FluxEPModelT0(lep)
            config!(epm_; verbose = true)
            converge!(epm_)
            return epm_
        end

        Q = MvNormal(epm)
        push!(Qs, Q)
    end

    # Original EP
    _, epm0 = withcachedat(PROJ, :get!, (:FluxEPModelT0, lep0,)) do 
        epm_ = FluxEPModelT0(lep0)
        converge!(epm_)
        return epm_
    end
    global Q0 = MultivariateNormal(epm0)

    return nothing
end

## ------------------------------------------------------------------
# MixtureModel

## ------------------------------------------------------------------
@time let
    
    global _unif_mix = MixtureModel(Qs)
    global _maxHu_mix = _max_mixture_H(
        _unif_mix, _entropy_ub1, GradientDescent()
    )
    # global _maxHl_mix = _max_mixture_H(
    #     _unif_mix, _entropy_lb1, GradientDescent()
    # )
end

## ------------------------------------------------------------------
let
    n = round(Int, 1e6)
    # @show _entropy_mcf(_unif_mix, n) # -15.290634234573337
    # @show _entropy_mcf(_maxHu_mix, n) # -1.58129779186625
    @show entropy(Q0)
    return nothing
end

## ------------------------------------------------------------------
let
    rxni = rand(1:35)
    p = plot(; xlabel = "flux [mmol/ gDW h]", ylabel = "pdf")
    _plot!(p, _unif_mix, rxni; label = "unif mix", lw = 3, c = :red)
    _plot!(p, _maxHu_mix, rxni; label = "maxHu mix", lw = 3, c = :blue)
    _plot!(p, Q0, rxni; label = "original", lw = 3, c = :black)
    # _plot!(p, _maxHl_mix, rxni; label = "maxHl mix", lw = 2)
    p
end

## ------------------------------------------------------------------
# marginal entropy
let
    idxs = 1:35
    _maxHu_marg_H = _entropy_mcf.(_marginal.([_maxHu_mix], idxs), Int(1e5));
    _unif_marg_H = _entropy_mcf.(_marginal.([_unif_mix], idxs), Int(1e5));
    _Q0_marg_H = _entropy_mcf.(_marginal.([Q0], idxs), Int(1e5));

    p = plot()
    plot!(p, _maxHu_marg_H; label = "maxHu", lw = 2)
    plot!(p, _unif_marg_H; label = "unif", lw = 2)
    plot!(p, _Q0_marg_H; label = "Q0", lw = 2)
    p
end
## ------------------------------------------------------------------
# MC
# let
#     global _mix = MixtureModel(Qs)

#     ns = 10.0.^[1:0.5:7;]
#     Hmcf = zeros(length(ns))
#     @threads :dynamic for (i, n) in collect(enumerate(ns))
#         @show n
#         Hmcf[i] = _entropy_mcf(_mix, round(Int, n))
#     end
#     plot(ns, Hmcf)
# end

## ------------------------------------------------------------------
function _inv(Σ)
    invΣ = similar(Σ)
    MetXBase.inplaceinverse!(invΣ, Σ)
    return invΣ
end

# ------------------------------------------------------------------
let
    # μ = (Σ₁⁻¹ + Σ₂⁻¹)⁻¹ (Σ₁⁻¹ μ₁ + Σ₂⁻¹ μ₂)
    # Σ = (Σ₁⁻¹ + Σ₂⁻¹)⁻¹
    invΣs = []
    μs = []
    
    for Q in Qs
        invΣ = _inv(Q.Σ)
        push!(invΣs, invΣ)
        push!(μs, Q.μ)
    end
    
    # Σ = (Σ₁⁻¹ + Σ₂⁻¹)⁻¹
    # sum(invΣs) = (Σ₁⁻¹ + Σ₂⁻¹)
    # Σ = _inv(sum(invΣs))
    Σ0 = _inv(sum(invΣs))
    μ0 = Σ0 * sum(invΣs .* μs)

    global Qprod = MvNormal(μ0, Σ0)
    nothing
end

## ------------------------------------------------------------------
function _plot_marginal!(p, Q0, i; 
        d = 3, kwargs...
    )
    μ = mean(Q0)[i]
    σ = sqrt(var(Q0)[i])

    N = Normal(μ, σ)
    TN = N
    # TN = Truncated(N, lb(epm, biom_id), ub(epm, biom_id))

    plot!(p, (x) -> pdf(TN, x), μ .- d*σ, μ .+ d*σ;
        kwargs...
    )
end

function _plot_marginal_sum!(p, Qs, i; d = 5, kwargs...)
    μs = [mean(Q)[i] for Q in Qs]
    σs = [sqrt(var(Q)[i]) for Q in Qs]

    xs = range(
        minimum(μs) .- d*maximum(σs),
        maximum(μs) .+ d*maximum(σs); 
        length = 1000
    )

    Qsum = zeros(length(xs))
    for (μ, σ) in zip(μs, σs)
        N = Normal(μ, σ)
        Qsum .+= pdf.([N], xs)
    end
    Z = sum(Qsum .* step(xs))
    plot!(p, xs, Qsum ./ Z; kwargs...)
end

## ------------------------------------------------------------------
let
    # lep0 EP
    # EP
    _, epm0 = withcachedat(PROJ, :get!, (:FluxEPModelT0, lep0,)) do 
        epm_ = FluxEPModelT0(lep0)
        converge!(epm_)
        return epm_
    end
    Q0 = MultivariateNormal(epm0)

    for ridx in 1:99
        p = plot(; xlabel = "flux", ylabel = "pdf")
        _plot_marginal!.([p], Qs, ridx; label = "", lw = 3, c = :red, alpha = 0.5)
        _plot_marginal!(p, Qs[1], ridx; label = "Q lim", lw = 3, c = :red, alpha = 0.5)
        # _plot_marginal!(p, Qprod, ridx; label = "prod Q lim", lw = 3, c = :orange)
        _plot_marginal_sum!(p, Qs, ridx; label = "sum Q lim", lw = 3, c = :blue)
        _plot_marginal!(p, Q0, ridx; label = "Q original", lw = 3, c = :black)
        
        sfig(PROJ, p, 
            ["HEK293"], 
            "Q_marginals", (;ridx, Qs = length(Qs)), ".png"
        ) |> println
    end


end

## ------------------------------------------------------------------


## ------------------------------------------------------------------
## ------------------------------------------------------------------
# @show mean(epm, biom_id)
# @show var(epm, biom_id)
# @show untruncated_mean(epm, biom_id)
# @show untruncated_var(epm, biom_id)
# N = Normal(
#     untruncated_mean(epm, biom_id), 
#     sqrt(untruncated_var(epm, biom_id))
# )
# TN = Truncated(N, lb(epm, biom_id), ub(epm, biom_id))

# plot((x) -> pdf(TN, x), 0.04, 0.05)

# N * N