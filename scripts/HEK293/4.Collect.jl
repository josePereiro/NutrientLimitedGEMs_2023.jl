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
    DB = TagDB()

    sims_dir = procdir(PROJ, ["HEK293", "sims"])

    sims_count = 0

    @time for (fni, fn) in enumerate(readdir(sims_dir; join = true))
        
        endswith(fn, ".jls") || continue
        _, sim = ldat(fn)
        sim["status"] == :success || continue
        haskey(sim, SIM_ID) || continue
        
        epdat = sim[SIM_ID]

        # stuff
        Fs = epdat["Fs"]
        Ss = epdat["Ss"]
        ep_statuses = epdat["ep_statuses"]

        # val_idxs1 = findall(ep_statuses .== :converged)
        val_idxs1 = last(ep_statuses) == :converged
        val_idxs2 = NutrientLimitedGEMs._find_val_idxs(Ss, Fs)
        val_idxs = intersect(val_idxs1, val_idxs2)
        isempty(val_idxs) && continue

        sims_count += 1

        # Collect
        simid = basename(fn)

        # @show fni, fn
    end

    @show sims_count

end


## ------------------------------------------------------------------
## ------------------------------------------------------------------
## ------------------------------------------------------------------
# frec
function __plot_ko_rxn_frec(PROJ::Project0, netid, ep_alg_version)
    
    net0 = pull_net("Martinez_Monge_HEK293");
    lep0 = nothing;
    perm = nothing
    ps = Plots.Plot[]

    # figid = "ΔF_high_biom"
    figid = "ΔF"
    
    vals_0 = nothing

    for (lth, uth) in [
            # # Biomass
            # (-Inf, Inf),
            # (0.3, 0.4),
            # (0.2, 0.3),
            # (0.1, 0.2),
            # (0.0, 0.1),

            # ΔF
            (-Inf, Inf),
            (-60.0, -50.0),
            (-50.0, -40.0),
            (-40.0, -30.0),
            (-30.0, -20.0),
            (-20.0, -10.0),
            (-10.0, -0.0),
        ]

        println("="^60)
        @show lth uth

        all_traj_idxs = []
        all_lp_sols = []
        ΔF1s = []

        _foreach_sim(PROJ, "HEK293", SIM_ID) do sim, epdat, val_idxs

            # biom
            biom1 = sim["biom1"]
            # push!(biom1s, biom1)
            # biom1 >= 0.2 || return nothing
            # biom1 <= 0.4 || return nothing
            
            Fs = epdat["Fs"]
            ΔFs = Fs[val_idxs] .- first(Fs[val_idxs])
            ΔFs *= -1
            last(ΔFs) >= lth || return nothing
            last(ΔFs) <= uth || return nothing
            push!(ΔF1s, last(ΔFs))

            # LP
            lep = sim["lep"]
            opm = FBAOpModel(lep, LP_SOLVER)
            optimize!(opm)
            sol = solution(opm)
            push!(all_lp_sols, sol)

            # traj_lens
            traj_idxs = sim["traj_idxs"]
            push!(all_traj_idxs, traj_idxs)

            if isnothing(lep0)
                lep0 = sim["lep0"]
            end
        
            return false

        end

        subSystems = subsystems(net0, colids(lep0))

        @show size(lep0)
        @show size(net0)

        C = length(all_traj_idxs)
        @show C
        
        C == 0 && continue

        # align_mat = zeros(Int, C, size(lep0, 2))
        # for (neti, idxs) in enumerate(all_traj_idxs)
        #     align_mat[neti, idxs] .= 1
        # end
        sols_mat = zeros(Float64, C, size(lep0, 2))
        for (neti, sol) in enumerate(all_lp_sols)
            sols_mat[neti, :] .= sol
        end

        # vals_ = count.(!iszero, eachcol(align_mat))
        vals_ = mean.(!iszero, eachcol(sols_mat))
        if isnothing(vals_0)
            vals_0 = vals_
            continue
        end
        vals_ = vals_ .- vals_0
        # if isnothing(perm)
            # perm = sortperm(vals_; rev = true)[1:90]
            perm = sortperm(abs.(vals_); rev = true)[1:90]
            perm = sort(perm; 
                rev = true,
                by = (i) -> begin
                    subSyst_map[subSystems[i]]
                end
            )
        # end
        
        # rxn_ids = colids(lep0, perm)
        # rxn_ids = [
        #     string("[", subSyst_map[subSystems[idx]], "] ", lep0.rxns[idx]) 
        #     for idx in perm
        # ]
        # rxn_ids = [subSyst_map[subSystems[idx]]
        #     for idx in perm
        # ]
        rxn_ids = [subSystems[idx]
            for idx in perm
        ]

        colors = [
            subSyst_colors[subSyst_map[subSystems[idx]]]
            for idx in perm
        ]

        xs = eachindex(perm)
        # ys = vals_[perm] ./ maximum(vals_[perm])
        ys = vals_[perm]
        p = bar(xs, ys;
            xrotation = 90,
            xticks = (xs, rxn_ids),
            size = (2*900, 500),
            yaxis = false,
            margin = 35Plots.mm,
            c = colors,
            bar_width = 1,
            label = "",
            title = string("filter: ", (lth, uth), " count: ", C),
            # tickfontsize = 10
            tickfontsize = 6
        )
        push!(ps, p)
        
    end
    
    sfig(PROJ, ps, 
        [netid], "ko_rxn_lp_sol", figid, ep_alg_version, ".png";
        layout = (length(ps), 1)
    ) |> println

end

# ------------------------------------------------------------------
let
    __plot_ko_rxn_frec(PROJ, "HEK293", SIM_ID) 
    # nothing
end

## ------------------------------------------------------------------
## ------------------------------------------------------------------
# Trajectoreis
function __plot_ko_rxn_heatmap(PROJ::Project0, netid, ep_alg_version)
    
    net0 = pull_net("Martinez_Monge_HEK293");
    lep0 = nothing;
    perm = nothing
    ps = Plots.Plot[]

    # figid = "ΔF_high_biom"
    figid = "ΔF"
    
    vals_0 = nothing
    
    all_traj_idxs = []
    all_lp_sols = []
    ΔF1s = []
    biom1s = []

    _foreach_sim(PROJ, "HEK293", SIM_ID) do sim, epdat, val_idxs

        # biom
        biom1 = sim["biom1"]
        push!(biom1s, biom1)

        Fs = epdat["Fs"]
        ΔFs = Fs[val_idxs] .- first(Fs[val_idxs])
        ΔFs *= -1
        push!(ΔF1s, last(ΔFs))

        # LP
        lep = sim["lep"]
        opm = FBAOpModel(lep, LP_SOLVER)
        optimize!(opm)
        sol = solution(opm)
        push!(all_lp_sols, sol)

        # traj_lens
        traj_idxs = sim["traj_idxs"]
        push!(all_traj_idxs, traj_idxs)

        if isnothing(lep0)
            lep0 = sim["lep0"]
        end
    
        return false

    end

    subSystems = subsystems(net0, colids(lep0))
    
    C = length(all_traj_idxs)
    @show C

    @show size(lep0)
    @show size(net0)

    # align_mat = zeros(Int, C, size(lep0, 2))
    # for (neti, idxs) in enumerate(all_traj_idxs)
    #     align_mat[neti, idxs] .= 1
    # end

    sols_mat = zeros(Float64, C, size(lep0, 2))
    for (neti, sol) in enumerate(all_lp_sols)
        sols_mat[neti, :] .= sol
    end

    return lep0, sols_mat

    # vals_ = count.(!iszero, eachcol(align_mat))
    vals_ = mean.(!iszero, eachcol(sols_mat))
    if isnothing(vals_0)
        vals_0 = vals_
        return
    end
    vals_ = vals_ .- vals_0
    # if isnothing(perm)
        # perm = sortperm(vals_; rev = true)[1:90]
        perm = sortperm(abs.(vals_); rev = true)[1:90]
        perm = sort(perm; 
            rev = true,
            by = (i) -> begin
                subSyst_map[subSystems[i]]
            end
        )
    # end
    # xs = eachindex(perm)
    # rxn_ids = colids(lep0, perm)
    # rxn_ids = [
    #     string("[", subSyst_map[subSystems[idx]], "] ", lep0.rxns[idx]) 
    #     for idx in perm
    # ]
    # rxn_ids = [subSyst_map[subSystems[idx]]
    #     for idx in perm
    # ]
    # rxn_ids = [subSystems[idx]
    #     for idx in perm
    # ]
    # colors = [
    #     subSyst_colors[subSyst_map[subSystems[idx]]]
    #     for idx in perm
    # ]

    sfig(PROJ, p, 
        [netid], "ko_rxn_lp_sol", figid, ep_alg_version, ".png";
        layout = (length(ps), 1)
    ) |> println

end

# ------------------------------------------------------------------
global lep0_, sols_mat_0 = __plot_ko_rxn_heatmap(PROJ, "HEK293", SIM_ID) 

## ------------------------------------------------------------------
let
    sols_mat_ = deepcopy(sols_mat_0)
    stds = sqrt.(var.(eachcol(sols_mat_)))
    cols_perm = sortperm(stds)

    # Normalize
    for coli in cols_perm
        sols_mat_[:, coli] .-= mean(sols_mat_[:, coli])
    end

    heatmap(view(sols_mat_, :, cols_perm))

    plot(log10.(stds[cols_perm] .+ 1e-6))
end