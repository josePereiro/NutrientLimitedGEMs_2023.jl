# This script is for creating the human networks from Chemostat_Rath2017.jl uuid = "f511193e-9926-11ea-07ad-fb434f320181"
@time begin

    using COBREXA
    using MAT
    using Serialization
    using ProjAssistant
    
    using NutrientLimitedGEMs
    using NutrientLimitedGEMs: load_HartGEM_model
    const NL = NutrientLimitedGEMs
    using GLPK
    using MetNets # this is not a dep

    using Base.Threads
end

## ------------------------------------------------------------------
# tools
function _make_unique(f::Function, vec)

    len = lastindex(vec)
    for (i, val) in enumerate(vec)

        for j in (i+1):len
            vali = vec[j]
            val == vali && f(j)
        end

    end

    return vec
end

_make_unique!(vec::Vector{String}) = _make_unique((i) -> vec[i] = string(vec[i], "_", i), vec)

## ------------------------------------------------------------------
# COBREXA hacks
"""
    metabolite_charge(m::MATModel, mid::String)::Maybe{Int}

Extract metabolite charge from `metCharge` or `metCharges`.
"""
function COBREXA.metabolite_charge(m::MATModel, mid::String)::Maybe{Int}
    met_charge = COBREXA._maybemap(
        x -> x[findfirst(==(mid), metabolites(m))],
        gets(m.mat, nothing, COBREXA._constants.keynames.metcharges),
    )
    isnothing(met_charge) ? 
        nothing :
        COBREXA._maybemap(Int, isnan(met_charge) ? nothing : met_charge)
end

## ------------------------------------------------------------------
function _to_cobrexa_CoreModel(srcfile)

    old_model = deserialize(srcfile)[:dat]
    
    M, N = size(old_model)
    _make_unique!(old_model.rxns)
    _make_unique!(old_model.mets)
    @assert allunique(old_model.rxns)
    @assert allunique(old_model.mets)
    return CoreModel(
        old_model.S,
        #= b =# zeros(M),
        #= c =# zeros(N),
        #= lb =# old_model.lb,
        old_model.ub,
        old_model.rxns,
        old_model.mets,
    )

end

# ------------------------------------------------------------------
_fill_if_illdim(x, dim, deflt) = length(x) == dim ? x : fill(deflt, dim)

function _to_COBREAXA_StandardModel(srcfile)

    met_readible = lrawdat(NL, "met_readable_ids.jls")
    exch_met_map = lrawdat(NL, "exch_met_map.jls")
    
    old_model = deserialize(srcfile)[:dat]
    M, N = size(old_model.S)
    _make_unique!(old_model.rxns)
    _make_unique!(old_model.mets)

    @assert allunique(old_model.rxns)
    @assert allunique(old_model.mets)

    # to dict
    dict_model = Dict()
    dict_model["S"] = old_model.S
    dict_model["c"] = _fill_if_illdim(Vector(old_model.c), N, 0.0)
    dict_model["lb"] = Vector(old_model.lb)
    dict_model["rxnNames"] = _fill_if_illdim(string.(old_model.rxnNames), N, "")
    dict_model["subSystems"] = _fill_if_illdim(string.(old_model.subSystems), N, "")
    dict_model["b"] = _fill_if_illdim(Vector(old_model.b), M, 0.0)
    
    # reformat
    dict_model["metCompartments"] = [
        string(met[end] == 's' ? "e" : met[end]) 
        for met in old_model.mets
    ]
    dict_model["mets"] = [
        string(met[1:end - 1], "[", dict_model["metCompartments"][i], "]") 
        for (i, met) in enumerate(old_model.mets)
    ]
    dict_model["metNames"] = [
        get(met_readible, met, "")
        for met in old_model.mets
    ]

    _format_exch = (id) -> begin
        id = string(id)
        haskey(exch_met_map, id) || return id
        startswith(id, "EX_") && return id
        return string("EX_", id)
    end
    dict_model["rxns"] = [
        _format_exch(rxn)
        for rxn in old_model.rxns
    ]

    dict_model["ub"] = Vector(old_model.ub)

    tempfile = replace(srcfile, ".jls" => ".mat")
    @info("Saving temp", tempfile)
    # show info
    for (k, val) in dict_model
        println(k, ": ", size(val))
    end
    
    matwrite(tempfile, Dict("model" => dict_model))

    # global model = load_model(file)
    @info("Re-Loading")
    model = load_model(tempfile)
    for (k, val) in model.mat
        eltype(val) == Any || continue
        model.mat[k] = string.(model.mat[k])
    end
    @info "Converting"
    model = convert(StandardModel, model)
    rm(tempfile; force = true)

    return model
end

## ------------------------------------------------------------------
# Process all .jls (They must be missing)
let
    head = "HartGEM"
    files = collect(readdir(rawdir(NL); join = true))
    for srcfile in files
        thid = threadid()

        endswith(srcfile, ".jls") || continue
        contains(srcfile, head) || continue

        for modelid in ["base", "fva_base", "fva_scaled", "scaled", "stst_scaled"], 
            tissue in ["GBM"], stst in ["A", "E"]

            if contains(srcfile, "=$modelid") && contains(srcfile, "=$tissue") && contains(srcfile, "=$stst")
                destfile = rawdir(NL, head, (;modelid, tissue, stst), ".json")
            elseif contains(srcfile, "=$modelid") && contains(srcfile, "=$tissue")
                destfile = rawdir(NL, head, (;modelid, tissue), ".json")
            else
                continue
            end
            
            # cache
            isfile(destfile) && continue

            println("\n")
            println("-"^60)
            @info("Doing", src=basename(srcfile), dest = basename(destfile), thid)
            model = _to_COBREAXA_StandardModel(srcfile)
            change_objective!(model, HUMAN_BIOMASS_IDER)

            # reformat exchanges
            
            # test FBA
            @info("Testing FBA", thid)
            coremodel = convert(CoreModel, model)
            sol = flux_balance_analysis_vec(coremodel, GLPK.Optimizer)
            biomidx = findfirst(HUMAN_BIOMASS_IDER .== coremodel.rxns)
            biomub = coremodel.xu[biomidx]
            biom = sol[biomidx]

            if biom > 0.01
                @info("FBA worked", biom, biomub, thid)
                @info("Saving", thid)
                save_model(model, destfile)
            else
                error("FBA failed!!! biom: ", biom, " ub: ", biomub)
            end

        end # for modelid 
    end # for srcfile
end

## ------------------------------------------------------------------
# Test loading
@time load_HartGEM_model(;modelid="fva_scaled", tissue="GBM")


