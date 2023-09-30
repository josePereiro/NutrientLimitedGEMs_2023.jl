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
# load ensembels
let
    _phase_ids = [
        "ECOLI-CORE-BEG2007-PHASE_1",
        "ECOLI-CORE-BEG2007-PHASE_2",
        "ECOLI-CORE-BEG2007-PHASE_3",
    ]
    _ensems_fns = Dict(
        "ECOLI-CORE-BEG2007-PHASE_1" => "12.1_ensem_ph1_zU.jl...<<len=503>>.jls",
        # "ECOLI-CORE-BEG2007-PHASE_2" => "",
        "ECOLI-CORE-BEG2007-PHASE_3" => "12.1_ensem_ph3_zU.jl...<<len=548>>.jls",
    )
    global ENSEMS_DAT = Dict()
    for (ph, fname) in _ensems_fns
        fn = procdir(PROJ, ["ensembles"], fname)
        _, dat = ldat(fn)
        ENSEMS_DAT[ph] = dat
    end
    nothing
end

## --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
function _microarray_patt(rxnid)
    # load micorarray data
    array_db = query(["ROOT", "PHASES_MICROARRAY"])
    global RAW_EXP_MAT = array_db["exp_mat"]
    haskey(RAW_EXP_MAT, rxnid) || nothing 
    # return time_tags = array_db["time_tags"]

    _patt = Float64
    # @show RAW_EXP_MAT[id]
    _mat = sum(RAW_EXP_MAT[rxnid]; dims = 1)
    # ph1: 2.0, 2.5, 3.5 
    # ph2: 4.0, 4.5, 5.0, 5.5, 6.0
    # ph3: 6.5, 7.0, 7.5, 8.0

    # else
    # ph_idxs = 
end

## --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# create ph expression and ense flux mats
let
    

    RAW_EXP_MAT["ACALD"]
    # for (rxn, genes) in raw_map["Supp3"]
    #     isempty(genes) && continue
    #     nrows = length(genes)
    #     ncols = length(genes[1]["raw"])
    #     exp_mat[rxn] = zeros(nrows, ncols)
    #     for (ri, obj) in enumerate(genes)
    #         exp_mat[rxn][ri, :] .= obj["raw"]
    #     end
    # end

    
end
