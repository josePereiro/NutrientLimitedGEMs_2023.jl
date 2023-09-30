# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
@time begin
    using Random
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
# TODO: create a IderVector, just a vector (Ordered Dict) which
# implements the ider interface
## --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# ensemble ph1 v1
# 1. Only glc allowed 
# 2. Average biomass fixed (sampled from an Uniform distribution)
let
    # context
    _simver = "ECOLI-CORE-BEG2007-PHASE_1"
    _load_contextdb(_simver)

    # lep
    CORE_XLEP_DB = query(["ROOT", "CORE_XLEP"])
    core_elep0 = CORE_XLEP_DB["core_elep0"][]
    global core_lep0 = lepmodel(core_elep0)
    RIDX = Dict(id => i for (i, id) in enumerate(colids(core_lep0)))
    core_elep0 = nothing
    
    tries_per_bb = 2000
    ensem_size = 500
    global ensem = []
    
    # biomass distribution
    ave_biom = 0.60 # from Beg2007 fig2.a
    B = Uniform(ave_biom - 0.1, ave_biom + 0.1)
    biom_dist_target = rand(B, floor(Int, ensem_size * 1.2)) 
    biom_th = 0.05
    _target_mean_biom = mean(biom_dist_target)
    @show mean(biom_dist_target)

    # find an essemble
    n0 = 0
    n1 = Inf
    nthrs = 10
    lk = ReentrantLock()
    perm = (fns) -> sort!(fns; by = (fn) -> rand()) # shuffle
    info_time = -1.0
    info_frec = 10 # second
    gc_time = -1.0
    gc_frec = 5 * 60 # second
    while length(ensem) < ensem_size
        _th_readdir(_simver; n0, n1, nthrs, perm) do bbi, bb
            # info
            lock(lk) do
                info_flag = time() > info_time
                if info_flag 
                    println("[", getpid(), ".", threadid(), "] ", 
                        "length(ensem) ", length(ensem), "\\", ensem_size, " ",
                        basename(rootdir(bb))
                    )
                    info_time = time() + info_frec
                end
                gc_flag = time() > gc_time
                if gc_flag
                    GC.gc()
                    gc_time = time() + gc_frec
                end
            end
            
            haskey(bb["meta"], "core_biomass_fba.ver") || return :ignore
            length(ensem) > ensem_size && return :break

            # frames
            feasets_frame = bb["core_feasets"]

            # peek
            tries = 0
            while true
                tries += 1
                tries > tries_per_bb && return :ignore

                # unpack
                feasets_blob0 = rand(feasets_frame)
                fealen, feaobj = rand(feasets_blob0)
                _core_sol = feaobj["core_biomass_fba.solution"]
                isempty(_core_sol) && continue

                # --------------------------------------------
                # boolean filters
                # --------------------------------------------
                # intake patterns
                # At phase 1, only glucose can be consumed
                good_pattern = true
                for exch in [
                        "EX_lac__D_e", "EX_malt_e",
                        "EX_gal_e", "EX_glyc_e"
                    ]
                    flx = _core_sol[RIDX[exch]]
                    abs(flx) < 1e-2 && continue
                    good_pattern = false
                    break
                end
                good_pattern || continue
                
                # --------------------------------------------
                # distribution filters
                # --------------------------------------------
                # biomass distribution
                ixd = RIDX["BIOMASS_Ecoli_core_w_GAM"]
                _core_biom = _core_sol[ixd]
                lock(lk) do
                    _thidx = findfirst(biom_dist_target) do z0
                        abs(_core_biom - z0) < biom_th
                    end
                    isnothing(_thidx) && return true
                    deleteat!(biom_dist_target, _thidx)
                    return false
                end && continue

                # --------------------------------------------
                # Push! obj
                # --------------------------------------------
                # Cool add ensemble
                lock(lk) do
                    push!(ensem, feaobj)
                end
                
            end # while true
        end # _th_readdir
        
        # increase tolerance
        biom_th = biom_th * 1.1
    end # while 

    # store
    _dat = Dict("ensem" => ensem, "RIDX" => RIDX)
    fn = procdir(PROJ, ["ensembles"], basename(@__FILE__), ".jls")
    sdat(_dat, fn; verbose = true)
    fn = procdir(PROJ, ["ensembles"], basename(@__FILE__), (;len = length(ensem)), ".jls")
    sdat(_dat, fn; verbose = true)

    println()
    println("= "^30)

    println("ENSEMBLE")
    println("- length(ensem)     ", length(ensem))
    
    flxs = _ensem_fba_solutions(core_lep0, ensem, "BIOMASS_Ecoli_core_w_GAM")
    println("- target mean(BIOM) ", _target_mean_biom)
    println("- ensem mean(BIOM)  ", mean(flxs))
    println("- ensem std(BIOM)   ", std(flxs))

    for each in [
            "EX_glc__D_e", "EX_lac__D_e", "EX_malt_e",
            "EX_gal_e", "EX_glyc_e", "EX_ac_e"
        ]
        flxs = _ensem_fba_solutions(core_lep0, ensem, each)
        println("- ensem mean($each)  ", mean(flxs))
        println("- ensem std($each)   ", std(flxs))
    end

    nothing
end

