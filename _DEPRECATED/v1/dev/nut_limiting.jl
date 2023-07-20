## ------------------------------------------------------------------
## ------------------------------------------------------------------
@time begin
    using MetXBase
    using MetXOptim
    using MetXNetHub
    using MetXEP
    using MetXCultureHub
    using Gurobi
    using Serialization
    
    using Plots
    using MetXPlots
    
    # Pkg.add("https://github.com/josePereiro/ImgTools.jl")
    using ImgTools
end

## ------------------------------------------------------------------
# TODO: Move to #TUTORIALS
# v^2 dual price
let

    global model_id = "ecoli_core"
    global net = pull_net(model_id)
    clampbounds!(net, -1000.0, 1000.0)
    
    glc_id = extras(net, "EX_GLC")
    glu_id = extras(net, "EX_GLU")
    nh4_id = extras(net, "EX_NH4")
    biom_id = extras(net, "BIOM")
    biom_idx = rxnindex(net, biom_id)

    # contextualization
    lb!(net, glc_id, -10.0)
    lb!(net, glu_id, -10.0)
    lb!(net, nh4_id, -10.0)

    # max biom
    global opm = FBAFluxOpModel(net, SOLVER)
    optimize!(opm)
    sol0 = solution(opm)
    biom0 = solution(opm, biom_id)

    # glc: double prices

    npoints = 5
    v0 = lb(net, glc_id)
    v1 = v0 * 0.9
    @assert v0 < 0 && v1 < 0
    test_points = range(v0, v1; length = npoints)

    # open competidors
    # lb!(net, glc_id, -10.0)
    lb!(net, glu_id, -100.0) 
    lb!(net, nh4_id, -100.0)
    
    # FBAFluxOpModel
    # bound biom
    ub!(net, biom_id, biom0)
    global opm = FBAFluxOpModel(net, SOLVER)
    ms0, errs = flux_dual_prices(opm, glc_id, test_points)
    glc_biom_m, glc_biom_err = ms0[biom_idx], errs[biom_idx]
    @show glc_biom_m, glc_biom_err
    
    p = plot(; title = model_id, xlabel = "flx idx", ylabel = "dual price")

    idxs = sortperm(abs.(ms0); rev = true)
    scatter!(p, ms0[idxs]; 
        # yerr = errs[idxs], 
        label = "(max z)", 
        # lw = 2, 
        c = :blue, 
        alpha = 0.6
    )

    # R2FBAFluxOpModel
    # fix biom
    bounds!(net, biom_id, biom0, biom0)
    global opm = R2FBAFluxOpModel(net, SOLVER)
    ms1, errs = flux_dual_prices(opm, glc_id, test_points)
    glc_biom_m, glc_biom_err = ms1[biom_idx], errs[biom_idx]
    @show glc_biom_m, glc_biom_err
    
    scatter!(p, ms1[idxs]; 
        label = "(min v'v)", 
        c = :red, 
        alpha = 0.6
    )

    _diff = ms0 .- ms1
    plot!(p, _diff[idxs]; 
        label = "diff", 
        lw = 2, c = :black, 
        ls = :dot
    )

    _max = max.(abs.(ms0), abs.(ms1))
    plot!(p, _max[idxs]; 
        label = "abs max", 
        lw = 2, c = :green, 
        ls = :dot
    )

    _miss_tol = 1e-5
    _miss = (abs.(ms0) .> _miss_tol .&& abs.(ms1) .< _miss_tol)
    plot!(p, _miss[idxs] .* _diff; 
        label = "miss", 
        lw = 2, c = :gold, 
    )

    @show net.rxns[_miss]
    
    p

end


## ------------------------------------------------------------------
## ------------------------------------------------------------------
## ------------------------------------------------------------------
let 

    global net0 = _setup_heknet()
    # global net = pull_net("ecoli_core")
    biom_id = extras(net0, "BIOM")
    exchs_ids0 = [
        "EX_GLC", "EX_NH4", "EX_TYR","EX_HIS",
        # "EX_ASN","EX_ARG",
        # "EX_ILE","EX_PHE","EX_MET","EX_THR","EX_LYS","EX_LEU",
        # "EX_VAL","EX_TRP","EX_GLU","EX_GLN"
    ]
    exchs_ids = extras.([net0], exchs_ids0)
    exchs_idxs = rxnindex.([net0], exchs_ids)
    
    p = plot()
    for exch_id in exchs_idxs

        
        bs = range(0.1, 1.0; length = 10)
        Ss = Float64[]
        for b in bs

            @time try
                lb!(net0, exch_id, -b)

                global net = box(net0, SOLVER; 
                    verbose = true, 
                    reduce = true,
                    protect_obj = true, 
                )

                println("\n", "="^30)
                @show exch_id
                @show b 

                opm = FBAFluxOpModel(net, SOLVER)
                optimize!(opm)
                biom1 = solution(opm, biom_id)
                @show biom1

                S = _entropy(net)
                @show S
                push!(Ss, S)

                
            catch err
                @show err
                push!(Ss, NaN)
            end
            
            lb!(net0, exch_id, -1.0)

        end

        plot!(p, bs, Ss; label = exch_id)
    end
    p

end

## ------------------------------------------------------------------








































## ------------------------------------------------------------------
## ------------------------------------------------------------------
## ------------------------------------------------------------------
# LIMITING NUTRIENT
let
    global model_id = "ecoli_core"
    global net = pull_net("ecoli_core")
    global citekey = "seniorRegulationNitrogenMetabolism1975"
    global cul = pull_cul(citekey)
    
    # contextualize
    glc_id = extras(net, "EX_GLC")
    glu_id = extras(net, "EX_GLU")
    nh4_id = extras(net, "EX_NH4")
    biom_id = extras(net, "BIOM")
    biom_idx = rxnindex(net, biom_id)

    # culid = "GLC_Limited_GLU"
    culid = "GLU_Limited"
    rep = 5
    D = queryfirst(cul, "D", culid, rep; extract = "val")
    X = queryfirst(cul, "X", culid, rep; extract = "val")
    c_glc = queryfirst(cul, "c_glc", culid, rep; extract = "val")
    c_glu = queryfirst(cul, "c_glu", culid, rep; extract = "val")
    @show c_glc, c_glu
    glc_L = -(c_glc * D / X)
    glu_L = -(c_glu * D / X)
    @show glc_L, glu_L
    
    # lb!(net, glc_id, glc_L)
    # lb!(net, glu_id, glu_L)
    # lb!(net, nh4_id, 0.0)

    lb!(net, glc_id, -50.0)
    lb!(net, glu_id, -10.0)
    lb!(net, nh4_id, -50.0)

    # dual price
    npoints = 5
    global opm = FBAFluxOpModel(net, Gurobi.Optimizer)

    # glc/biom double price
    v0 = lb(net, glc_id)
    v1 = v0 * 0.9
    @assert v0 < 0 && v1 < 0

    test_points = range(v0, v1; length = npoints)
    ms, errs = lb_dual_prices(opm, glc_id, test_points)
    glc_biom_m, glc_biom_err = ms[biom_idx], errs[biom_idx]
    @show glc_biom_m, glc_biom_err

    # glu/biom double price
    v0 = glu_L
    v1 = v0 * 0.9
    @assert v0 < 0 && v1 < 0

    test_points = range(v0, v1; length = npoints)
    ms, errs = lb_dual_prices(opm, glu_id, test_points)
    glu_biom_m, glu_biom_err = ms[biom_idx], errs[biom_idx]
    @show glu_biom_m, glu_biom_err

    if abs(glu_biom_m) > abs(glc_biom_m)
        @info "GLU is limiting"
    else
        @info "GLC is limiting"
    end

    return nothing


end