## ------------------------------------------------------------------
begin
    using MetXNetHub
    using MetXOptim
    using MetXCultureHub
    using ContextDBs
    using Gurobi
end

## ------------------------------------------------------------------
function _test_fba(net)
    # biomass_human
    opm = FBAOpModel(net, Gurobi.Optimizer)
    optimize!(opm)
    @show solution(opm, "biomass_human")
    return 
end

## ------------------------------------------------------------------
# TODO: Think about a global id map package
function _exch_met_map()
    id_map = Dict()
    id_map["glc"] = "HMR_9034"
    id_map["gln"] = "HMR_9063"
    id_map["val"] = "HMR_9046"
    id_map["cys"] = "HMR_9065"
    id_map["pro"] = "HMR_9068"
    id_map["asn"] = "HMR_9062"
    id_map["ile"] = "HMR_9039"
    id_map["leu"] = "HMR_9040"
    id_map["arg"] = "HMR_9066"
    id_map["his"] = "HMR_9038"
    id_map["lys"] = "HMR_9041"
    id_map["met"] = "HMR_9042"
    id_map["phe"] = "HMR_9043"
    id_map["trp"] = "HMR_9045"
    id_map["tyr"] = "HMR_9064"
    id_map["ala"] = "HMR_9061"
    id_map["gly"] = "HMR_9067"
    id_map["ser"] = "HMR_9069"
    id_map["thr"] = "HMR_9044"
    id_map["glu"] = "HMR_9071"
    id_map["asp"] = "HMR_9070"
    id_map["amm"] = "EX_nh4[e]"
    id_map["galc"] = "HMR_9140"
    id_map["pyr"] = "HMR_9133"
    id_map["lac"] = "HMR_9135"

    for (k, v) in id_map
        id_map[v] = k
    end

    return id_map

end


## ------------------------------------------------------------------
# Set Rath medium bounds
function _set_rath_medium_bounds!(net, cul_db, culid)

    # culture data
    stst_db = queryall(cul_db, ["ROOT", "culid" => culid])
    c_db = queryall(stst_db, ["ROOT", "apiid" => r"c_"])
    @assert all(c_db[:, "unit"] .== "mM")
    D_db = query(stst_db, ["ROOT", "apiid" => "D"])
    D = D_db["val"]
    Xv_db = query(stst_db, ["ROOT", "apiid" => "Xv"])
    Xv = Xv_db["val"]

    # close exchanges
    for ri in eachindex(net.rxns)
        net.subSystems[ri] == "Exchange/demand reactions" || continue
        bounds!(net, ri, 0.0, 1000.0)
    end
    
    # Open medium
    # Medium from Rath 
    
    id_map = _exch_met_map()
    for rath_id in c_db[!, "metid"]
        # intake = D * c / Xv
        c = query(c_db, ["ROOT", "metid" => rath_id])["val"]
        l = -D * c / Xv
        # @show u
        
        gem_id = id_map[rath_id]
        bounds!(net, gem_id, l, 1000.0) 
    end

    # completed with the Ham medium from
    # https://github.com/SysBioChalmers/Human-GEM/data/metabolicTasks/metabolicTasks_Essential.txt GR Growth on Ham's media (biomass production) task

    bounds!(net, "HMR_9145", -1000.0, 1000.0) # pantothenate (metid: MAM02680e)
    bounds!(net, "HMR_9167", -1000.0, 1000.0) # lipoic acid (metid: MAM02394e)
    bounds!(net, "HMR_9074", -1000.0, 1000.0) # sulfate (metid: MAM02946e)
    bounds!(net, "HMR_9404", -1000.0, 1000.0) # retinoate (metid: MAM02833e)
    bounds!(net, "HMR_9035", -1000.0, 1000.0) # linoleate (metid: MAM02387e)
    bounds!(net, "HMR_9036", -1000.0, 1000.0) # linolenate (metid: MAM02389e)
    bounds!(net, "HMR_9083", -1000.0, 1000.0) # choline (metid: MAM01513e)
    # bounds!(net, "HMR_9361", -1000.0, 1000.0) # inositol (metid: MAM02171e)
    bounds!(net, "HMR_9143", -1000.0, 1000.0) # riboflavin (metid: MAM02842e)
    bounds!(net, "HMR_9378", -1000.0, 1000.0) # nicotinamide (metid: MAM02583e)
    bounds!(net, "HMR_9144", -1000.0, 1000.0) # pyridoxine (metid: MAM02817e)
    # bounds!(net, "HMR_9423", -1000.0, 1000.0) # thymidine (metid: MAM02996e)
    bounds!(net, "HMR_9109", -1000.0, 1000.0) # biotin (metid: MAM01401e) 
    bounds!(net, "HMR_9151", -1000.0, 1000.0) # alpha-tocopherol (metid: MAM01327e)
    bounds!(net, "HMR_9146", -1000.0, 1000.0) # folate (metid: MAM01830e)
    bounds!(net, "HMR_9358", -1000.0, 1000.0) # hypoxanthine (metid: MAM02159e)
    bounds!(net, "HMR_9269", -1000.0, 1000.0) # aquacob(III)alamin (metid: MAM01361e)
    bounds!(net, "HMR_9153", -1000.0, 1000.0) # gamma-tocopherol (metid: MAM01935e)
    bounds!(net, "HMR_9159", -1000.0, 1000.0) # thiamin (metid: MAM02982e)
    bounds!(net, "HMR_9048", -1000.0, 1000.0) # O2 (metid: MAM02630e)
    bounds!(net, "HMR_9047", -1000.0, 1000.0) # H2O (metid: MAM02040e)
    bounds!(net, "HMR_9072", -1000.0, 1000.0) # Pi (metid: MAM02751e)
    bounds!(net, "HMR_9076", -1000.0, 1000.0) # Fe2+ (metid: MAM01821e)

   nothing

end
