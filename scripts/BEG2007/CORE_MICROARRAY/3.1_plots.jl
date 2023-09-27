## ------------------------------------------------------------
@time begin
    using MetXGEMs
    using MetXBase
    using MetXNetHub
    using MetXOptim
    using ProjFlows
    using BlobBatches
    using Statistics
    using CairoMakie
    using NutrientLimitedGEMs
end

# TODO: Add Graphs.jl kind of functionality for getting basic stuff
# ------------------------------------------------------------
include("1_setup.jl")
include("2_utils.jl")

## ------------------------------------------------------------
# All phases microarray expression patterns
let
    # context
    array_db = query(["ROOT", "PHASES_MICROARRAY"])
    exp_mat = array_db["exp_mat"]
    time_tags = array_db["time_tags"]
    
    for (rxn, mat) in exp_mat
        f = Figure()
        ax = Axis(f[1,1]; 
            xlabel = "time [h]", ylabel = "expression level", 
            title = rxn
        )
        for patt in eachrow(mat)
            lines!(ax, time_tags, patt)
        end
        sfig(PROJ, f, [SIMVER, "microarrays"], 
            rxn, ".png"
        ) |> println
    end
end

## ------------------------------------------------------------
using Clp
let
    # context
    array_db = query(["ROOT", "PHASES_MICROARRAY"])
    global exp_mat = array_db["exp_mat"]
    global time_tags = array_db["time_tags"]

    # net
    global net0 = pull_net("ecoli_core_Beg2007")
    global net = box(net0, Clp.Optimizer)
end

## ------------------------------------------------------------
let
    opm = fba(net0, Clp.Optimizer)

    corcoes = Float64[]
    for timei in 1:12
        flxs = Float64[]
        rexps = Float64[]
        for (rxn, mat) in exp_mat
            if !hascolid(net0, rxn) 
                @show rxn
                continue
            end
            for patt in eachrow(mat)
                flx = solution(opm, rxn)
                rexp = patt[timei] ./ maximum(patt)
                push!(flxs, abs(flx))
                push!(rexps, rexp)
            end
        end
        push!(corcoes, cor(flxs, rexps))
    end

    f = Figure()
    # ax = Axis(f[1,1]; xlabel = "abs fluxes", ylabel = "rel expression")
    ax = Axis(f[1,1]; xlabel = "time [h]", ylabel = "correlation", title = "original model")
    barplot!(ax, time_tags, corcoes)
    f
end

## ------------------------------------------------------------
# Comparizon with downregulated ensembles
let

    # context
    _load_contextdb(SIMVER)
    array_db = query(["ROOT", "PHASES_MICROARRAY"])
    exp_mat = array_db["exp_mat"]
    time_tags = array_db["time_tags"]

    phase = "ECOLI-CORE-BEG2007-PHASE_I-0.1.0"
    _load_contextdb(phase)

     # context
    GLOB_DB = query(["ROOT", "GLOBALS"])
    XLEP_DB = query(["ROOT", "CORE_XLEP"])
    LP_SOLVER = GLOB_DB["LP_SOLVER"]
    DOWNREG_FACTOR = GLOB_DB["DOWNREG_FACTOR"] # TOSYNC

    global f = Figure()
    ax = Axis(f[1,1])

    global _max_m = 0.0

    _th_readdir(phase; n0 = 1, n1 = 3, nthrs = 1) do bbi, bb

        # top filter
        haskey(bb["meta"], "core_nut_sp.ver") || return :ignore

        # lep
        core_elep0 = XLEP_DB["core_elep0"][]
        core_lep0 = lepmodel(core_elep0)
        core_elep0 = nothing

        # opm
        opm = FBAOpModel(core_lep0, LP_SOLVER)

        # frames
        feasets_frame = bb["core_feasets"]
        strip_frame = bb["core_strip"]

        # run
        info_frec = 50
        for (blobi, (strip_blob, feasets_blob)) in enumerate(
                zip(strip_frame, feasets_frame)
            )
            
            # info
            info_flag = blobi == 1 || blobi == lastindex(feasets_frame) || iszero(rem(blobi, info_frec)) 
            info_flag && println("[", getpid(), ".", threadid(), "] ", 
                "blobi ", blobi, "\\", length(feasets_frame), " ",
                basename(rootdir(bb))
            )

            # TODO: test correlations with downregulation vectors against delta expression
            
            koset = strip_blob["koset"]
            for (li, feaobj) in feasets_blob
                
                feaset = koset[1:li]
                biom_max = feaobj["core_biomass_fba.biom"]
                fba_sol = feaobj["core_biomass_fba.solution"]
                glc_sp = get(feaobj, "core_nut_sp.EX_glc__D_e.obj_m", 0.0)
                lac_sp = get(feaobj, "core_nut_sp.EX_lac__D_e.obj_m", 0.0) # max 0.042877197265625
                malt_sp = get(feaobj, "core_nut_sp.EX_malt_e.obj_m", 0.0) # max 0.17822265625
                gal_sp = get(feaobj, "core_nut_sp.EX_gal_e.obj_m", 0.0)
                glyc_sp = get(feaobj, "core_nut_sp.EX_glyc_e.obj_m", 0.0)
                ac_sp = get(feaobj, "core_nut_sp.EX_ac_e.obj_m", 0.0)

                _max_m = max(_max_m, abs(lac_sp))

                # filter
                # biom_max < 0.7 || continue
                # biom_max > 0.5 || continue
                # abs(glc_sp) > 0.04 || continue
                # abs(lac_sp) > 0.04 || continue
                # abs(malt_sp) > 0.04 || continue
                # abs(gal_sp) > 0.04 || continue
                # abs(glyc_sp) < 1e-3 || continue

                # correlation with microarray data
                corcoes = Float64[]
                for timei in 1:12
                    flxs = Float64[]
                    rexps = Float64[]
                    for (rxn, mat) in exp_mat
                        hascolid(core_lep0, rxn) || continue
                        rxni = colindex(core_lep0, rxn)
                        patt = sum.(eachcol(mat))
                        # for patt in eachrow(mat)
                            
                            # flx = fba_sol[rxni]
                            flx = ub(core_lep0, rxni) - lb(core_lep0, rxni)
                            push!(flxs, abs(flx))
                            # rflx = flx > 0 ? 
                            #     flx / ub(core_lep0, rxni) : 
                            #     flx / lb(core_lep0, rxni) 
                            # rflx = isnan(rflx) ? 0.0 : rflx
                            # push!(flxs, abs(rflx))
                            
                            # rexp = patt[timei] ./ maximum(patt)
                            rexp = patt[timei]
                            push!(rexps, rexp)
                        # end
                    end
                    push!(corcoes, cor(flxs, rexps))
                end

                # plot
                lines!(ax, time_tags, corcoes)

            end # for li
            
        end # for blobi
    end
    f
end