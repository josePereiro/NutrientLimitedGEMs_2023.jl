# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
@time begin
    using Random
    using MetXGEMs
    using MetXBase
    using Statistics
    using MetXEP
    using BlobBatches
    using CairoMakie
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
# TODO: create a IderVector, just a vector (Ordered Dict) which
# implements the ider interface
## --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# biomass ensemble
let
    # context
    _simver = "ECOLI-CORE-BEG2007-PHASE_I-0.1.0"
    _load_contextdb(_simver)

    # lep
    CORE_XLEP_DB = query(["ROOT", "CORE_XLEP"])
    core_elep0 = CORE_XLEP_DB["core_elep0"][]
    global core_lep0 = lepmodel(core_elep0)
    core_elep0 = nothing
    
    ensem_size = 100
    global ensem = []
    tries_per_bb = 100
    # biomass
    ave_biom = 0.6
    ave_biom1 = 1.2
    biom_th = 1e-2
    B = _MaxEnt_beta(ave_biom)
    biom_dist_target = rand(B, ensem_size) .* ave_biom1
    @show mean(biom_dist_target)

    # find an essemble
    n0 = 0
    n1 = Inf
    nthrs = 1
    while !isempty(biom_dist_target)
        perm = (fns) -> sort!(fns; by = (fn) -> rand()) # shuffle
        _th_readdir(;n0, n1, nthrs, perm) do bbi, bb
            @show bb
            @show length(ensem)
            
            haskey(bb["meta"], "core_biomass_fba.ver") || return :ignore

            # frames
            feasets_frame = bb["core_feasets"]
            # sort!(feasets_frame; by = (fn) -> rand()) # shuffle

            # peek
            tries = 0
            while true
                tries += 1
                tries > tries_per_bb && return :ignore

                # unpack
                feasets_blob0 = rand(feasets_frame)
                fealen, feaobj = rand(feasets_blob0)
                _core_sol = feaobj["core_biomass_fba.solution"]

                # --------------------------------------------
                # boolean filters
                # --------------------------------------------
                # biomass range
                
                # intake patterns
                # At phase 1, only glucose can be consumed
                good_pattern = true
                for exch in [
                        "EX_lac__D_e", "EX_malt_e",
                        "EX_gal_e", "EX_glyc_e"
                    ]
                    ixd = colindex(core_lep0, exch)
                    abs(_core_sol[ixd] < 1e-2) || continue
                    good_pattern = false
                end
                good_pattern || continue

                
                # --------------------------------------------
                # distribution filters
                # --------------------------------------------
                # biomass distribution
                ixd = colindex(core_lep0, "BIOMASS_Ecoli_core_w_GAM")
                _core_biom = _core_sol[ixd]
                _idx = findfirst(biom_dist_target) do z0
                    abs(_core_sol - z0) < biom_th
                end
                isnothing(_idx) && continue

                # TODO: Store necessary data about the net
                push!(ensem, feaobj)
                deleteat!(biom_dist_target, _idx)
                isempty(biom_dist_target) && return :break
            end
        end # _th_readdir
        
        # increase tolerance
        biom_th = biom_th * 1.1
    end # while 

    # store
    

    nothing
end

