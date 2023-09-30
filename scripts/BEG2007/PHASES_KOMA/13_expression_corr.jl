# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
@time begin
    using MetXGEMs
    using MetXBase
    using Statistics
    using MetXEP
    using BlobBatches
    using CairoMakie
    using Distributions
    using MetXOptim
    using Statistics
    using MetXEP
    using Clp
    using Base.Threads
    using NutrientLimitedGEMs
end

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
include("1_setup.jl")
include("2_utils.jl")

## --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
function _microarray_patts()
    # load micorarray data
    array_db = query(["ROOT", "PHASES_MICROARRAY"])
    RAW_EXP_MAT = array_db["exp_mat"]
    MICRO_PATTS = Dict()
    for rxnid in keys(RAW_EXP_MAT)
        _patt = sum(RAW_EXP_MAT[rxnid]; dims = 1)
        _patt .= _patt ./ mean(_patt)
        _patt = [mean(_patt[1:4]), mean(_patt[5:8]), mean(_patt[9:12])]
        MICRO_PATTS[rxnid] = _patt
    end
    return MICRO_PATTS
end

function _ensems_mean_vec(fnhint::Regex)
    dir = procdir(PROJ, ["ensembles"])
    files = readdir(dir; join = true)
    sort!(files; by = f -> rand()) # shuffle

    ENS_MEANS = Dict{String, Float64}()
    for fn in files
        m = match(fnhint, basename(fn))
        isnothing(m) && continue
        
        _, dat = ldat(fn)
        ensem = dat["ensem"]
        for (rxn, idx) in dat["RIDX"]
            flxs = _ensem_fba_solutions(ensem, idx)
            ENS_MEANS[rxn] = mean(flxs)
        end
        break
    end
    return ENS_MEANS
end

function _ensems_patt(hits...)
    _means = map(_ensems_mean_vec, hits)
    # all_keys = string.(unique(vcat(collect.(keys.(_means))...)))
    _comm = intersect([Set(keys(_m)) for _m in _means]...)
    ENS_PATTS = Dict()
    for _m in _means
        for key in _comm
            _patt = get!(ENS_PATTS, key, Float64[])
            push!(_patt, _m[key])
        end
    end
    ENS_PATTS
end
_ensems_patt(r"PHASE_1.*Uniform.*5000")
# ECOLI-CORE-BEG2007-PHASE_1...Uniform...5000...<<h=10001632393434744400>>.jls
## --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# create ph expression and ense flux mats
let
    

    
end
