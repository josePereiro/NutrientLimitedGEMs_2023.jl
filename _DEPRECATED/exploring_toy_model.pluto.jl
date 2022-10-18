### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 010bfa4a-033f-11ed-3204-013a734525b6
# ╠═╡ show_logs = false
let
	import Pkg
	proj = Base.current_project(@__DIR__)
	Pkg.activate(proj)
end

# ╔═╡ 9e185c8e-53bf-4bf2-91e3-c2cc97d23c9c
# ╠═╡ show_logs = false
begin
    using Plots
	using PlutoUI
    using COBREXA
    using NutrientLimitedGEMs
    const NL = NutrientLimitedGEMs
    using GLPK
    using ProjAssistant
end

# ╔═╡ ec018419-8777-4e00-ad19-72dbfe7c5782
begin
	model0 = NL.load_toy_model3();
	abs_max = 100
	bound_range = 0.0:0.01:abs_max
	u_bounds = zeros(length(model0.xu))
	l_bounds = zeros(length(model0.xu))
end;

# ╔═╡ db73f85d-bd91-49f2-8554-a819cc5df47d
# modify bounds
let
	# "Ex_s", "Ex_a", "Ex_l", "vspe", "vap", "vas", "vpe", "vpl", "biomass"
	# model.xu
	# model.xl
	# model.rxns
end

# ╔═╡ fd066a64-e6b8-432e-aa79-3e535b916041
let
	
	md"""

	### Model bounds
	
	Ex [s->] ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` `
	lb:  $(@bind Ex_s_lb Slider(bound_range; default = abs(model0.xl[1])))
	` ` ` ` ` `
	ub:  $(@bind Ex_s_ub Slider(bound_range; default = model0.xu[1])) \\\
	Ex [a->] ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` `
	lb:  $(@bind Ex_a_lb Slider(bound_range; default = abs(model0.xl[2])))
	` ` ` ` ` `
	ub:  $(@bind Ex_a_ub Slider(bound_range; default = model0.xu[2])) \\\
	Ex [l->] ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` 
	lb:  $(@bind Ex_l_lb Slider(bound_range; default = abs(model0.xl[3])))
	` ` ` ` ` `
	ub:  $(@bind Ex_l_ub Slider(bound_range; default = model0.xu[3])) \\\
	v [s->p+e] ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` 
	lb:  $(@bind vspe_lb Slider(bound_range; default = abs(model0.xl[4])))
	` ` ` ` ` `
	ub:  $(@bind vspe_ub Slider(bound_range; default = model0.xu[4])) \\\
	v [a->p]  ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` `
	lb:  $(@bind vap_lb Slider(bound_range; default = abs(model0.xl[5])))
	` ` ` ` ` `
	ub:  $(@bind vap_ub Slider(bound_range; default = model0.xu[5])) \\\
	v [a->s]  ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` `
	lb:  $(@bind vas_lb Slider(bound_range; default = abs(model0.xl[6])))
	` ` ` ` ` `
	ub:  $(@bind vas_ub Slider(bound_range; default = model0.xu[6])) \\\
	v [p->e]  ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` `
	lb:  $(@bind vpe_lb Slider(bound_range; default = abs(model0.xl[7])))
	` ` ` ` ` `
	ub:  $(@bind vpe_ub Slider(bound_range; default = model0.xu[7])) \\\
	v [p->l]  ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` ` `
	lb:  $(@bind vpl_lb Slider(bound_range; default = abs(model0.xl[8])))
	` ` ` ` ` `
	ub:  $(@bind vpl_ub Slider(bound_range; default = model0.xu[8])) \\\
	biomass [a + e + s ->] ` ` ` ` ` ` ` ` ` ` ` `
	lb:  $(@bind z_lb Slider(bound_range; default = abs(model0.xl[9])))
	` ` ` ` ` `
	ub:  $(@bind z_ub Slider(bound_range; default = model0.xu[9])) \\\
		
	"""
end

# ╔═╡ f7c5414f-a157-4a8a-be66-bfb76bd135bf
begin
	model = deepcopy(model0)
	model.xl .= -[Ex_s_lb, Ex_a_lb, Ex_l_lb, vspe_lb, vap_lb, vas_lb, vpe_lb, vpl_lb, z_lb]
	model.xu .= [Ex_s_ub, Ex_a_ub, Ex_l_ub, vspe_ub, vap_ub, vas_ub, vpe_ub, vpl_ub, z_ub]
	@assert all(model.xu .>= model.xl)
	model
end;

# ╔═╡ 2bc9f6c7-890a-44a4-8f67-474f7b7c85b1


# ╔═╡ e417d1a4-4616-4f23-9cc4-7f13a8028aad
begin
	dat = Dict{String, Any}()
end;

# ╔═╡ eba1aa65-b11e-41bf-a686-1e07fa7b87ff
let
	model

	ps = Plots.Plot[]
	
	p = plot()
	scatter!(p, reverse(model.xu), reverse(model.rxns); label = "", ms = 5, c = :red)
	scatter!(p, reverse(model.xl), reverse(model.rxns); label = "", ms = 5, c = :blue)
	plot!(p; xlabel = "bounds", size = [200, 400])
	push!(ps, p)

	for exchid in find_uptake_exchs(model) 

		samples, mat = dat[exchid]
		exchidx = rxnidx(model, exchid)
		vals = mat[exchidx, :]
		
		p = plot()
		plot!(p, vals, mat[1, :]; label = "Ex_s", lw = 3)
		plot!(p, vals, mat[2, :]; label = "Ex_a", lw = 3)
		plot!(p, vals, mat[end, :]; label = "biomass", lw = 3, ls = :dash)
		plot!(;xlabel = exchid, ylabel = "flux", legend = false)
		push!(ps, p)
	end

	layout = @layout [
    	a{0.3w} [b
        	c]
	]
	plot(ps...; size = [800, 500], margin = 3Plots.mm, layout)
end

# ╔═╡ 2d73206d-b6d2-4ad1-82f6-d1404e5d9e31


# ╔═╡ e9a0b7d8-1532-49e7-88ac-f5001036c14a
let
	nsample = 100
	# compute dep mat
    for exchid in find_uptake_exchs(model)

        # samples
        lb, ub = bounds(model, exchid)
        lb, ub = max(lb, -abs_max), min(abs(lb), abs_max)
        samples = lb .+ (range(0.0, 1.0; length = nsample) .* (ub - lb))
        
        # compute
        vals, mat = compute_effect_mat(
            model, exchid, samples;
			clearcache = true, 
			verbose = false
        )
		dat[exchid] = (vals, mat)
    end
end

# ╔═╡ ae977152-25dc-42e0-a843-455993648893
# ider  $(@bind ider Select(model.rxns, default = "biomass"))
# 	beta  $(@bind beta Slider(D["beta_range"]))
# 	alpha $(@bind alpha Slider(D["alpha_range"]; 
# 		default = maximum(D["alpha_range"])))


# ╔═╡ 39cdcea1-6bd5-4cc5-a1eb-dfca908c7656


# ╔═╡ Cell order:
# ╟─010bfa4a-033f-11ed-3204-013a734525b6
# ╟─9e185c8e-53bf-4bf2-91e3-c2cc97d23c9c
# ╟─ec018419-8777-4e00-ad19-72dbfe7c5782
# ╟─db73f85d-bd91-49f2-8554-a819cc5df47d
# ╟─fd066a64-e6b8-432e-aa79-3e535b916041
# ╟─f7c5414f-a157-4a8a-be66-bfb76bd135bf
# ╠═eba1aa65-b11e-41bf-a686-1e07fa7b87ff
# ╠═2bc9f6c7-890a-44a4-8f67-474f7b7c85b1
# ╟─e417d1a4-4616-4f23-9cc4-7f13a8028aad
# ╠═2d73206d-b6d2-4ad1-82f6-d1404e5d9e31
# ╠═e9a0b7d8-1532-49e7-88ac-f5001036c14a
# ╠═ae977152-25dc-42e0-a843-455993648893
# ╠═39cdcea1-6bd5-4cc5-a1eb-dfca908c7656
