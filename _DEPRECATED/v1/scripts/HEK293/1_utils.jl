## ------------------------------------------------------------------
begin
    using NutrientLimitedGEMs_2023
    using ProjFlows
    using MetX
    using MetXBase
    using Gurobi
    using DataFrames
    using Plots
    using Random
    using Distributions
    using Base.Threads
    using Optim
    using LinearAlgebra
end

## ------------------------------------------------------------------
# Experimental data
# TODO: package Calzadilla_et_al
_, (CalDat, _) = lprocdat(PROJ,
    ["Calzadilla_et_al"], "Calzadilla.bundle", ".jls"; 
    verbose = false
)
nothing

## ------------------------------------------------------------------
_Calzadilla_hek_idmap = Dict(
    "ALA" => "R_EX_ala_L_LPAREN_e_RPAREN_",
    "ARG" => "R_EX_arg_L_LPAREN_e_RPAREN_",
    "ASN" => "R_EX_asn_L_LPAREN_e_RPAREN_",
    "ASP" => "R_EX_asp_L_LPAREN_e_RPAREN_",
    "GLN" => "R_EX_gln_L_LPAREN_e_RPAREN_",
    "GLU" => "R_EX_glu_L_LPAREN_e_RPAREN_",
    "GLC" => "R_EX_glc_LPAREN_e_RPAREN_",
    "GLY" => "R_EX_gly_LPAREN_e_RPAREN_",
    "ILE" => "R_EX_ile_L_LPAREN_e_RPAREN_",
    "LAC" => "R_EX_lac_L_LPAREN_e_RPAREN_",
    "LEU" => "R_EX_leu_L_LPAREN_e_RPAREN_",
    "LYS" => "R_EX_lys_L_LPAREN_e_RPAREN_",
    "MET" => "R_EX_met_L_LPAREN_e_RPAREN_",
    "PHE" => "R_EX_phe_L_LPAREN_e_RPAREN_",
    "PRO" => "R_EX_pro_L_LPAREN_e_RPAREN_",
    "SER" => "R_EX_ser_L_LPAREN_e_RPAREN_",
    "THR" => "R_EX_thr_L_LPAREN_e_RPAREN_",
    "TRP" => "R_EX_trp_L_LPAREN_e_RPAREN_",
    "TYR" => "R_EX_tyr_L_LPAREN_e_RPAREN_",
    "VAL" => "R_EX_val_L_LPAREN_e_RPAREN_",
)
for (k,v) in _Calzadilla_hek_idmap
    _Calzadilla_hek_idmap[v] = k
end

## ------------------------------------------------------------------
function _Calzadilla_heknet!(lep0, st;
        solver = LP_SOLVER,
    )

    # contextualization
    # lb = - c D / Xv

    stdat = CalDat[st]
    D = stdat["D"]["val"]
    Xv = stdat["Xv"]["val"]

    # ["ALA", "ARG", "ASN", "ASP", "GLN", "GLU", "GLC", 
    # "GLY", "ILE", "LAC", "LEU", "LYS", "MET", "PHE", 
    # "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
    for CalId in CalDat["metids"]

        c = stdat["c$(CalId)"]["val"]
        lb = -(c * D / Xv)
        
        HekId = _Calzadilla_hek_idmap[CalId]
        # println(CalId, " lb: ", lb)
        lb!(lep0, HekId, lb)

    end
    
    # TODO: ask Calzadilla
    lb!(lep0, "R_EX_tyr_L_LPAREN_e_RPAREN_", -1.0) # To restrictive
    lb!(lep0, "R_EX_arg_L_LPAREN_e_RPAREN_", -1.0) # To restrictive
    lb!(lep0, "R_EX_his_L_LPAREN_e_RPAREN_", -1.0)
    
    lb!(lep0, "R_EX_o2_LPAREN_e_RPAREN_", -1000.0)
    lb!(lep0, "R_EX_h_LPAREN_e_RPAREN_", -1000.0)
    lb!(lep0, "R_EX_pi_LPAREN_e_RPAREN_", -1000.0)
    lb!(lep0, "R_EX_h2o_LPAREN_e_RPAREN_", -1000.0)
    lb!(lep0, "R_EX_co2_LPAREN_e_RPAREN_", -1000.0)
    lb!(lep0, "R_EX_hco3_LPAREN_e_RPAREN_", -1000.0)

    # max biom
    # opm = FBAOpModel(lep0, solver)
    # optimize!(opm)
    # biom_id = extras(lep0, "BIOM")
    # biom0 = solution(opm, biom_id)
    # @show biom0
    # @show D

    # _, lep = withcachedat(PROJ, :get!, (lep0,);
    #     verbose = true
    # ) do 
    #     box(lep0, solver; verbose = true)
    # end
    return lep0
end

function _Calzadilla_heknet(arg...; kwargs...)
    model_id = "Martinez_Monge_HEK293"
    net = pull_net(model_id)
    clampbounds!(net, -1000.0, 1000.0)
    _Calzadilla_heknet!(net, arg...; kwargs...)
end

## ------------------------------------------------------------------
# Martinez_Monge_HEK293
function _base_heknet(PROJ::Project0; 
        solver = LP_SOLVER,
    )

    model_id = "Martinez_Monge_HEK293"
    net = pull_net(model_id)
    clampbounds!(net, -1000.0, 1000.0)

    # old contextualization
    lb!(net, "R_EX_tyr_L_LPAREN_e_RPAREN_", -1.0)
    lb!(net, "R_EX_his_L_LPAREN_e_RPAREN_", -1.0)
    lb!(net, "R_EX_asn_L_LPAREN_e_RPAREN_", -1.0)
    lb!(net, "R_EX_arg_L_LPAREN_e_RPAREN_", -1.0)
    lb!(net, "R_EX_ile_L_LPAREN_e_RPAREN_", -1.0)
    lb!(net, "R_EX_phe_L_LPAREN_e_RPAREN_", -1.0)
    lb!(net, "R_EX_met_L_LPAREN_e_RPAREN_", -1.0)
    lb!(net, "R_EX_thr_L_LPAREN_e_RPAREN_", -1.0)
    lb!(net, "R_EX_lys_L_LPAREN_e_RPAREN_", -1.0)
    lb!(net, "R_EX_leu_L_LPAREN_e_RPAREN_", -1.0)
    lb!(net, "R_EX_val_L_LPAREN_e_RPAREN_", -1.0)
    lb!(net, "R_EX_trp_L_LPAREN_e_RPAREN_", -1.0)

    lb!(net, "R_EX_glc_LPAREN_e_RPAREN_", -1.0)

    lb!(net, "R_EX_glu_L_LPAREN_e_RPAREN_", -1.0)
    lb!(net, "R_EX_gln_L_LPAREN_e_RPAREN_", -1.0)
    # lb!(net, "R_EX_lac_L_LPAREN_e_RPAREN_", -1000.0)
    # lb!(net, "R_EX_nh4_LPAREN_e_RPAREN_", -1000.0)

    lb!(net, "R_EX_o2_LPAREN_e_RPAREN_", -1000.0)
    lb!(net, "R_EX_h_LPAREN_e_RPAREN_", -1000.0)
    lb!(net, "R_EX_pi_LPAREN_e_RPAREN_", -1000.0)
    lb!(net, "R_EX_h2o_LPAREN_e_RPAREN_", -1000.0)
    lb!(net, "R_EX_co2_LPAREN_e_RPAREN_", -1000.0)
    lb!(net, "R_EX_hco3_LPAREN_e_RPAREN_", -1000.0)

    # max biom
    # opm = FBAOpModel(net, solver)
    # optimize!(opm)
    # biom_id = extras(net, "BIOM")
    # biom0 = solution(opm, biom_id)
    # @show biom0
    # # @show D

    _, lep = withcachedat(PROJ, :get!, (net,);
        verbose = true
    ) do 
        box(net, solver; verbose = true)
    end
    return lep
end

## ------------------------------------------------------------------
function _foreach_sim(f::Function, PROJ, netid, simid)
    
    traj_dir = procdir(PROJ, [netid, "sims"])

    @time for fn in readdir(traj_dir; join = true)
        endswith(fn, ".jls") || continue
        _, sim = ldat(fn)
        sim["status"] == :success || continue
        haskey(sim, simid) || continue
        
        epdat = sim[simid]

        # stuff
        Fs = epdat["Fs"]
        Ss = epdat["Ss"]
        ep_statuses = epdat["ep_statuses"]

        val_idxs1 = findall(ep_statuses .== :converged)
        val_idxs2 = NutrientLimitedGEMs_2023._find_val_idxs(Ss, Fs)
        val_idxs = intersect(val_idxs1, val_idxs2)
        isempty(val_idxs) && continue

        f(sim, epdat, val_idxs) === true && break
    end
end

## ------------------------------------------------------------------
subSyst_map = Dict(
    "Transport, mitochondrial" => "Transport",
    "Glycine, serine, alanine and threonine metabolism" => "AA metabolism",
    "Lysine metabolism" => "AA metabolism",
    "Tryptophan metabolism" => "AA metabolism",
    "Sphingolipid metabolism" => "Lipid metabolism",
    "Valine, leucine, and isoleucine metabolism" => "AA metabolism",
    "Transport, extracellular" => "Transport",
    "Fatty acid oxidation" => "Lipid metabolism",
    "Citric acid cycle" => "Citric acid cycle",
    "Nucleotide interconversion" => "Nucleotide metabolism",
    "Purine synthesis" => "AA metabolism",
    "Triacylglycerol synthesis" => "Lipid metabolism",
    "Methionine and cysteine metabolism" => "AA metabolism",
    "null" => "",
    "Glutamate metabolism" => "AA metabolism",
    "Urea cycle" => "Urea cycle",
    "R group synthesis" => "AA metabolism",
    "Pyrimidine synthesis" => "Pyrimidine synthesis",
    "Alanine and aspartate metabolism" => "AA metabolism",
    "Oxidative phosphorylation" => "Oxidative phosphorylation",
    "Cholesterol metabolism" => "Lipid metabolism",
    "Glycerophospholipid metabolism" => "Lipid metabolism",
    "Transport, endoplasmic reticular" => "Transport",
    "Transport, nuclear" => "Transport",
    "Folate metabolism" => "Folate metabolism",
    "Exchange/demand reaction" => "Transport",
    "Squalene and cholesterol synthesis" => "Lipid metabolism",
    "Glycolysis/gluconeogenesis" => "Sugar metabolism",
    "" => "",
    "Fatty acid synthesis" => "Lipid metabolism",
    "Arginine and Proline Metabolism" => "AA metabolism",
    "Pentose phosphate pathway" => "Sugar metabolism",
    "Nucleotide salvage pathway" => "Nucleotide metabolism",
    "Histidine metabolism" => "AA metabolism",
    "Miscellaneous" => "Miscellaneous",
    "Pyruvate metabolism" => "Pyruvate metabolism",
    "Inositol phosphate metabolism" => "Inositol phosphate metabolism",
    "Purine catabolism" => "Nucleotide metabolism",
    "Cysteine Metabolism" => "AA metabolism",
)

subSyst_colors = begin 
    subs = unique(collect(values(subSyst_map)))
    colors = Plots.distinguishable_colors(length(subs))
    Dict(subs .=> colors)
end


# ------------------------------------------------------------------
# ------------------------------------------------------------------
# MixtureModel

## ------------------------------------------------------------------
function _rand_mixture(M::Int, D::Int = 1; 
        ws = ones(M), 
        μscale = 1.0, μ0 = 0.0,
        σscale = 1.0, σ0 = 0.0,
    )

    qs = MvNormal[]
    for _ in 1:M
        μ = μ0 .+ μscale .* randn(D)
        Σ = σ0 .+ σscale .* rand(D, D)
        Σ = Σ'Σ
        push!(qs, MvNormal(μ, Σ))
    end

    ws = abs.(ws)
    ws ./= sum(ws)

    return MixtureModel(qs, ws)
end

# ------------------------------------------------------------------
# Ploting tools
function _suggest_range(q::Normal, s)
    return (q.μ - s * q.σ, q.μ + s * q.σ)
end

function _suggest_range(q::MvNormal, s)
    return (q.μ .- s .* diag(q.Σ), q.μ .+ s .* diag(q.Σ))
end

function _suggest_range(nm::MixtureModel, s)
    x0, x1 = Inf, -Inf
    for q in nm.components
        _x0, _x1 = _suggest_range(q, s)
        x0 = min.(x0, _x0)
        x1 = max.(x1, _x1)
    end
    return x0, x1
end

function _marginal(mv::MvNormal, k::Int)
    return Normal(mean(mv)[k], sqrt(var(mv)[k]))
end

function _marginal(mm::MixtureModel, k::Int)
    _mqs = Normal[_marginal(q, k) for q in components(mm)]
    return MixtureModel(_mqs, probs(mm))
end

function _plot!(p::Plots.Plot, q::UnivariateDistribution, x0, x1; 
        nsamples = 1000, norm = identity,
        plot_kwargs...
    )
    xs = range(x0, x1; length = nsamples)
    ys = [pdf(q, x) for x in xs]
    plot!(p, xs, norm(ys); plot_kwargs...)
end

function _plot!(p::Plots.Plot, q::UnivariateDistribution; plot_kwargs...)
    x0, x1 = _suggest_range(q, 4)
    _plot!(p, q, x0, x1; plot_kwargs...)
end

_plot!(p::Plots.Plot, q::MultivariateDistribution, k::Int, args...; kwargs...) = 
    _plot!(p, _marginal(q, k), args...; kwargs...)

_plot(args...; kwargs...) = _plot!(plot(), args...; kwargs...)

# ------------------------------------------------------------------
# Entropy

# ------------------------------------------------------------------
# Bounds

function _z_factor(di::Normal, dj::Normal)
    nz = Normal(dj.μ, di.σ + dj.σ)
    return pdf(nz, di.μ)
end

function _z_factor(di::MvNormal, dj::MvNormal)
    nz = MvNormal(dj.μ, di.Σ + dj.Σ)
    return pdf(nz, di.μ)
end

function _entropy_lb1(mm::MixtureModel)
    
    # Hl = - \sum_i w_i \log( \sum_j w_j z_ij )
    # z_ij = N( μ_i | μ_j, C_i + C_j )

    M = ncomponents(mm)
    mix_qs = components(mm)
    mix_ws = probs(mm)

    # compute zij
    z = zeros(M, M)
    for i in 1:M, j in 1:M
        z[i, j] = _z_factor(mix_qs[i], mix_qs[j])
    end

    Hl = 0.0
    for i in 1:M
        sum_j = 0
        for j in 1:M
            sum_j += mix_ws[j] * z[i, j]
        end
        Hl += mix_ws[i] * sum_j
    end

    return -Hl
end

function _entropy_ub1(mm::MixtureModel)

    # Hu = \sum_i w_i ( - \log w_i + H(q))

    M = ncomponents(mm)
    mix_qs = components(mm)
    mix_ws = probs(mm)

    Hu = 0.0
    for i in 1:M
        Hu += mix_ws[i] * ( -log(mix_ws[i]) + entropy(mix_qs[i]))
    end
    return Hu
end

## ------------------------------------------------------------------
# simple quadrature
function _entropy_quad(
        mm::MixtureModel{<:UnivariateDistribution}, 
        ngrid = Int(1e5)
    )

    # H = \sum_{X} p(x) \log p(x)

    x0, x1 = _suggest_range(mm, 5)
    _xs = range(x0, x1; length = ngrid)
    _dx = step(_xs)
    _sum = 0.0
    for x in _xs
        p = pdf(mm, x)
        _sum += p * log(p)
    end
    return -1 * _sum * _dx
end

## ------------------------------------------------------------------
# Monte carlos
using Random
function _monte_carlo_integration_lazy(f, a, b, n)
    dx = prod(b .- a)
    # s = zeros(nthreads())
    s = 0.0
    box = Uniform.(a, b)
    # @threads :static for _ in 1:n
    for _ in 1:n
        s[threadid()] += f(rand.(box))
    end
    # return (1/n) * sum(s) * dx
    return (1/n) * s * dx
end

function _monte_carlo_integration_lazy(f, dist::Distribution, n::Int)
    # integral = zeros(nthreads())
    integral = 0.0
    # @threads :static for i in 1:n
    for i in 1:n
        sample = rand(dist)
        weight = pdf(dist, sample)
        # integral[threadid()] += f(sample) / weight
        integral += f(sample) / weight
    end
    # return sum(integral) / n
    return integral / n
end

function _entropy_mcu(mm::Distribution, a, b, n)
    H = _monte_carlo_integration_lazy(a, b, n) do x
        p = pdf(mm, x)
        return p * log(p)
    end
    return -H
end
function _entropy_mcu(mm::Distribution, n; s = 5)
    a, b = _suggest_range(mm, s)
    return _entropy_mcu(mm, a, b, n)
end

function _entropy_mcf(mm::Distribution, f::Distribution, n)
    H = _monte_carlo_integration_lazy(f, n) do x
        p = pdf(mm, x)
        return p * log(p)
    end
    return -H
end
_entropy_mcf(mm::Distribution, n) = _entropy_mcf(mm, mm, n)

## ------------------------------------------------------------------
# W Entropy maximization
function _max_mixture_H!(
        mm::MixtureModel, 
        Hfun::Function, 
        solver::Optim.AbstractOptimizer; 
        show_trace = true, 
        iterations = 1_000, 
        epsconv = 1e-30
    )

    mix_ws = probs(mm)
    
    f = (_ws) -> let 
        _ws .= abs.(_ws)
        _ws ./= sum(_ws)
        mix_ws .= _ws # update MixtureModel
        return -1 * Hfun(mm)
    end

    # Set the initial guess
    mix_ws .= abs.(mix_ws)
    mix_ws ./= sum(mix_ws)
    
    # Minimize the function
    opt = Optim.Options(; 
        show_trace, iterations, 
        g_tol = epsconv
    )
    result = optimize(f, mix_ws, solver, opt)
    # maxH_ws = result.minimizer
    # maxH = result.minimum

    mix_ws .= abs.(result.minimizer)
    mix_ws ./= sum(mix_ws)
    
    return mm
end

_max_mixture_H(mm::MixtureModel, Hfun::Function, solver::Optim.AbstractOptimizer; kwargs...) = 
    _max_mixture_H!(deepcopy(mm), Hfun, solver; kwargs...)

## ------------------------------------------------------------------
# Utils

## ------------------------------------------------------------------
function _rescale(v::Vector)
    _v = v .- mean(v)
    _v ./= sqrt(var(_v))
    return _v
end


nothing