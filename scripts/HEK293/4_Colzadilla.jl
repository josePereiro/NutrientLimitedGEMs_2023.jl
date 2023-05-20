## ------------------------------------------------------------------
@time begin
    using NutrientLimitedGEMs
    const NL = NutrientLimitedGEMs

    using ProjFlows
    using MetX
    using Gurobi
    using DataFrames
    using Plots
end

# ------------------------------------------------------------------
include("0_params.jl")
include("1_utils.jl")

## ------------------------------------------------------------------
let 
    global metids = CalDat["metids"]
    global met_colors_vec = Plots.distinguishable_colors(length(metids))
    global met_colors_dict = Dict((metids .=> met_colors_vec)...)
end

## ------------------------------------------------------------------
# Dev
# let
#     st = 4
    
#     stdat = CalDat[st]
#     D = stdat["D"]["val"]

#     lep1 = _Calzadilla_heknet(st)
#     biom_id = extras(lep1, "BIOM")

#     # test LP
#     println("> Calzadilla")
#     opm = FBAOpModel(lep1, LP_SOLVER)
#     optimize!(opm)
#     biom0 = solution(opm, biom_id)
#     @show biom0
#     @show D
#     println()
# end

## ------------------------------------------------------------------
# Loas nutlim nets
let
    global nl_leps = []
    traj_dir = procdir(PROJ, ["HEK293", "sims"])
    @time for fn in readdir(traj_dir; join = true)
        endswith(fn, ".jls") || continue
        _, sim = ldat(fn)
        sim["status"] == :success || continue
        haskey(sim, SIM_ID) || continue
        @show fn
        
        # box vol
        biom1 = sim["biom1"]
        biom1 > 0.1 || continue
        @show biom1
        
        lep = sim["lep"]        
        push!(nl_leps, lep)
        # return lep
    end
end


## ------------------------------------------------------------------
# Collect
global DB = TagDB()

# ------------------------------------------------------------------
@time let
    solver = LP_SOLVER
    
    # fixing_set = ["GLN", "LAC"]
    fixing_set = []
    fixing_sets = []
    for idxs in Iterators.product(fill([true, false], length(fixing_set))...)
        _subset = fixing_set[findall(idxs)]
        push!(fixing_sets, _subset)
    end
    
    for fixxing in fixing_sets

        for st in 1:5
            
            orig_net = _base_heknet(PROJ)
            
            biom_id = extras(orig_net, "BIOM")

            nets_pools = [
                ("original", [orig_net]), 
                ("nut_lim", nl_leps[1:25])
            ]

            for (netid0, nets_) in nets_pools
                
                for (neti, lep_) in enumerate(nets_)
                    
                    netid = string(netid0, neti)
                    lep0 = deepcopy(lep_)
                    
                    stdat = CalDat[st]
                    D = stdat["D"]["val"]
                    
                    println("-"^60)
                    @show st
                    @show netid
                    @show D
                    println()

                    try

                        # test LP
                        println("> Original")
                        opm = FBAOpModel(lep0, solver)
                        optimize!(opm)
                        biom0 = solution(opm, biom_id)
                        @show biom0
                        @show D
                        println()

                        # Contextualize
                        lep1 = _Calzadilla_heknet!(lep0, st)

                        # test LP
                        println("> Calzadilla")
                        opm = FBAOpModel(lep1, solver)
                        optimize!(opm)
                        biom0 = solution(opm, biom_id)
                        @show biom0
                        @show D
                        println()

                        # Fixing GLC and D
                        HEKId = _Calzadilla_hek_idmap["GLC"]
                        val = stdat["qGLC"]["val"]
                        dval = abs(val * 0.1)
                        bounds!(lep1, HEKId, val .- dval, val .+ dval)
                        
                        dval = abs(D * 0.2)
                        bounds!(lep1, biom_id, D .- dval, D .+ dval)

                        # Extra Fixing
                        for CalId in fixxing
                            HEKId = _Calzadilla_hek_idmap[CalId]
                            val = stdat["q$CalId"]["val"]
                            dval = abs(val * 0.1)
                            bounds!(lep1, HEKId, val .- dval, val .+ dval)
                        end

                        # test LP
                        println("> Fix D and GLC")
                        opm = FBAOpModel(lep1, solver)
                        optimize!(opm)
                        biom0 = solution(opm, biom_id)
                        @show biom0
                        @show D
                        println()

                        # boxing
                        _, lep = withcachedat(PROJ, :get!, (:Boxing, lep1,)) do 
                            box(lep1, solver; verbose = true)
                        end

                        # LP
                        println("> Boxed")
                        opm = FBAOpModel(lep, solver)
                        optimize!(opm)
                        biom0 = solution(opm, biom_id)
                        @show biom0
                        @show D
                        println()

                        # EP
                        _, epm = withcachedat(PROJ, :get!, (:FluxEPModelT0, lep,)) do 
                            epm_ = FluxEPModelT0(lep)
                            converge!(epm_)
                            return epm_
                        end

                        # Collect
                        push!(DB, 
                            "netid" => netid,
                            "st" => st, 
                            "fixxing" => fixxing,
                            "lep" => lep,
                            "epm" => epm,
                            "opm" => opm,
                        )
                        
                    catch err
                        # rethrow(err)
                        @warn err
                    end

                end # for st
            end # for nets
        end # for lep0
    end # fixxing

end

## ------------------------------------------------------------------
## ------------------------------------------------------------------
_mse(dat, ref) = median((dat .- ref).^2)

## ------------------------------------------------------------------
# i) Fix D, and GLC: Focus only into the big fluxes...
# ii) Explore the space of NutLim GEMs
# iii) Compute more trajectories.

let
    # netid0 = "original"
    netid0 = "nut_lim"
    objs = query(DB, 
        "netid" => Regex(string(netid0, ".*")), 
        # "st" => 1,
    )

    exp_vals = []
    exp_errs = []
    lp_vals = []
    ep_vals = []
    colors = []

    for obj in objs

        st = obj["st"]
        stdat = CalDat[st]
        fixxing = obj["fixxing"]
        length(fixxing) == 0 || continue
        # fixxing == ["GLC"] || continue
        
        # for CalId in ["GLN", "LAC", "VAL", "GLU"]
        # for CalId in ["LAC"]
        for CalId in stdat["metids"]
            HekId = _Calzadilla_hek_idmap[CalId]

            # err = get!(mse, CalId, Dict())
            # lperr = get!(err, CalId, Dict())
            # eperr = get!(err, CalId, Dict())
            
            exp_val = stdat["q$CalId"]["val"]
            exp_err = stdat["q$CalId"]["err"]
            
            lp_val = solution(obj["opm"], HekId)
            ep_val = mean(obj["epm"], HekId)
            
            push!(exp_vals, exp_val)
            push!(exp_errs, exp_err)
            push!(lp_vals, lp_val)
            push!(ep_vals, ep_val)
            push!(colors, met_colors_dict[CalId])
        end
    end

    sidxs = sortperm(abs.(exp_vals))

    # xs = eachindex(ep_vals)
    
    # p = plot(;
    #     xlabel = "dat index", 
    #     ylabel = "log flux", 
    #     title = netid0
    # )
    
    # scatter!(p, xs, abs.(ep_vals)[sidxs]; 
    #     label = "EP vs EXP", 
    #     m = (5, :square), 
    #     c = colors[sidxs]
    # )

    # scatter!(p, xs, abs.(lp_vals)[sidxs]; 
    #     label = "LP vs EXP", 
    #     m = (5, :circle), 
    #     c = colors[sidxs]
    # )

    # plot!(p, xs, 
    #     abs.(exp_vals)[sidxs]; 
    #     label = "EXP", 
    #     lw = 3, ls = :dash
    # )

    p = plot()
    scatter!(p,  exp_vals, lp_vals)
    scatter!(p, exp_vals, ep_vals)
end

## ------------------------------------------------------------------
