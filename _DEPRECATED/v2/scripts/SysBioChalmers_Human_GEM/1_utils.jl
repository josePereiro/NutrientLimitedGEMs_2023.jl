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
    @show solution(opm, "MAR13082")
    return 
end

## ------------------------------------------------------------------
# TODO: Think about a global id map package
function _exch_met_map()
    id_map = Dict()
    id_map["glc"] = "MAR09034"
    id_map["gln"] = "MAR09063"
    id_map["val"] = "MAR09046"
    id_map["cys"] = "MAR09065"
    id_map["pro"] = "MAR09068"
    id_map["asn"] = "MAR09062"
    id_map["ile"] = "MAR09039"
    id_map["leu"] = "MAR09040"
    id_map["arg"] = "MAR09066"
    id_map["his"] = "MAR09038"
    id_map["lys"] = "MAR09041"
    id_map["met"] = "MAR09042"
    id_map["phe"] = "MAR09043"
    id_map["trp"] = "MAR09045"
    id_map["tyr"] = "MAR09064"
    id_map["ala"] = "MAR09061"
    id_map["gly"] = "MAR09067"
    id_map["ser"] = "MAR09069"
    id_map["thr"] = "MAR09044"
    id_map["glu"] = "MAR09071"
    id_map["asp"] = "MAR09070"
    id_map["amm"] = "MAR11420"
    id_map["galc"] = "MAR09140"
    id_map["pyr"] = "MAR09133"
    id_map["lac"] = "MAR09135"

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

    bounds!(net, "MAR09145", -1000.0, 1000.0) # pantothenate (metid: MAM02680e)
    bounds!(net, "MAR09167", -1000.0, 1000.0) # lipoic acid (metid: MAM02394e)
    bounds!(net, "MAR09074", -1000.0, 1000.0) # sulfate (metid: MAM02946e)
    bounds!(net, "MAR09404", -1000.0, 1000.0) # retinoate (metid: MAM02833e)
    bounds!(net, "MAR09035", -1000.0, 1000.0) # linoleate (metid: MAM02387e)
    bounds!(net, "MAR09036", -1000.0, 1000.0) # linolenate (metid: MAM02389e)
    bounds!(net, "MAR09083", -1000.0, 1000.0) # choline (metid: MAM01513e)
    # bounds!(net, "MAR09361", -1000.0, 1000.0) # inositol (metid: MAM02171e)
    bounds!(net, "MAR09143", -1000.0, 1000.0) # riboflavin (metid: MAM02842e)
    bounds!(net, "MAR09378", -1000.0, 1000.0) # nicotinamide (metid: MAM02583e)
    bounds!(net, "MAR09144", -1000.0, 1000.0) # pyridoxine (metid: MAM02817e)
    # bounds!(net, "MAR09423", -1000.0, 1000.0) # thymidine (metid: MAM02996e)
    bounds!(net, "MAR09109", -1000.0, 1000.0) # biotin (metid: MAM01401e) 
    bounds!(net, "MAR09151", -1000.0, 1000.0) # alpha-tocopherol (metid: MAM01327e)
    bounds!(net, "MAR09146", -1000.0, 1000.0) # folate (metid: MAM01830e)
    bounds!(net, "MAR09358", -1000.0, 1000.0) # hypoxanthine (metid: MAM02159e)
    bounds!(net, "MAR09269", -1000.0, 1000.0) # aquacob(III)alamin (metid: MAM01361e)
    bounds!(net, "MAR09153", -1000.0, 1000.0) # gamma-tocopherol (metid: MAM01935e)
    bounds!(net, "MAR09159", -1000.0, 1000.0) # thiamin (metid: MAM02982e)
    bounds!(net, "MAR09048", -1000.0, 1000.0) # O2 (metid: MAM02630e)
    bounds!(net, "MAR09047", -1000.0, 1000.0) # H2O (metid: MAM02040e)
    bounds!(net, "MAR09072", -1000.0, 1000.0) # Pi (metid: MAM02751e)
    bounds!(net, "MAR09076", -1000.0, 1000.0) # Fe2+ (metid: MAM01821e)

   nothing

end
