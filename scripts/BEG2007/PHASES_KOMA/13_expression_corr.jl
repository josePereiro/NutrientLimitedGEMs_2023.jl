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
function _microarray_patts(;phs = 1:3)
    # load micorarray data
    array_db = query(["ROOT", "PHASES_MICROARRAY"])
    RAW_EXP_MAT = array_db["exp_mat"]
    MICRO_PATTS = Dict()
    for rxnid in keys(RAW_EXP_MAT)
        _patt = sum(RAW_EXP_MAT[rxnid]; dims = 1)
        _patt = [mean(_patt[1:4]), mean(_patt[5:8]), mean(_patt[9:12])]
        MICRO_PATTS[rxnid] = _patt[phs]
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
# ECOLI-CORE-BEG2007-PHASE_1...Uniform...5000...<<h=10001632393434744400>>.jls
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# create ph expression and ense flux mats
let
    _vs = Float64[]
    _Es = Float64[]
    _corrs = Float64[]
    for it in 1:10
        _ensems = _ensems_patt(r"PHASE_1.*Uniform.*5000", r"PHASE_3.*Uniform.*5000")
        _microarray = _microarray_patts(;phs = [1,3])
        _comm = intersect(Set(keys(_ensems)), Set(keys(_microarray)))
        __vs = Float64[]
        __Es = Float64[]
        for rxn in _comm
            for (v, E) in zip(_ensems[rxn], _microarray[rxn])
                push!(__vs, v)
                push!(__Es, E)
            end
        end
        push!(_corrs, cor(__vs, __Es))
        push!(_vs, __vs...)
        push!(_Es, __Es...)
    end

    f = Figure()
    ax = Axis(f[1,1];)
    _vs = abs.(_vs) .- mean(abs.(_vs))
    _vs ./= std(_vs)
    _Es = _Es .- mean(_Es)
    _Es ./= std(_Es)
    scatter!(ax, _vs, _Es; color = :black)
    lines!(ax, _vs, _vs; color = :black)
    @show cor(_vs, _Es)
    f
end
