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
    using NutrientLimitedGEMs_2023
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
# 2. Biomass \in [z0, z1]
let

    # lep
    CORE_XLEP_DB = query(["ROOT", "CORE_XLEP"])
    core_elep0 = CORE_XLEP_DB["core_elep0"][]
    core_lep0 = lepmodel(core_elep0)
    RIDX = Dict(id => i for (i, id) in enumerate(colids(core_lep0)))
    core_elep0 = nothing

    # local ARGS
    tries_per_bb = 500
    ENSEM_SIZE = parseARGS("ENS-SIZE:") do _str
        isnothing(_str) && return 5000
        return parse(Int, _str)
    end
    BIOMAS_RANGE = parseARGS("BIOMAS-RANGE:") do _str
        isnothing(_str) && return 0.0:1000.0
        return eval(Meta.parse(_str))::AbstractRange
    end

    return BIOMAS_RANGE
    
    while true

        ensem = []

        # find an essemble
        n0 = 0
        n1 = Inf
        nthrs = 5
        lk = ReentrantLock()
        perm = (fns) -> sort!(fns; by = (fn) -> rand()) # shuffle
        info_time = -1.0
        info_frec = 10 # second
        gc_time = -1.0
        gc_frec = 5 * 60 # second
        while length(ensem) < ENSEM_SIZE
            _th_readdir(SIMVER; n0, n1, nthrs, perm) do bbi, bb
                # info
                lock(lk) do
                    info_flag = time() > info_time
                    if info_flag 
                        println("[", getpid(), ".", threadid(), "] ", 
                            "length(ensem) ", length(ensem), "\\", ENSEM_SIZE, " ",
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
                length(ensem) > ENSEM_SIZE && return :break

                # frames
                feasets_frame = bb["core_feasets"]

                # peek
                tries = 0
                while true
                    tries += 1
                    tries > tries_per_bb && return :ignore

                    # unpack
                    feasets_blob0 = rand(feasets_frame)
                    isempty(feasets_blob0) && continue
                    fealen, feaobj = rand(feasets_blob0)
                    _core_sol = feaobj["core_biomass_fba.solution"]
                    isempty(_core_sol) && continue

                    # --------------------------------------------
                    # boolean filters
                    # --------------------------------------------
                    if SIMVER == "ECOLI-CORE-BEG2007-PHASE_0"
                        # do nothing
                    elseif SIMVER == "ECOLI-CORE-BEG2007-PHASE_1"
                        # do nothing
                    elseif SIMVER == "ECOLI-CORE-BEG2007-PHASE_2"
                        # intake patterns
                        # At phase 2 no EX_glyc_e nor EX_ac_e intake
                        good_pattern = true
                        for exch in ["EX_glyc_e"]
                            flx = _core_sol[RIDX[exch]]
                            if flx < -1e-2 
                                good_pattern = false
                                break
                            end
                        end
                        good_pattern || continue
                    elseif SIMVER == "ECOLI-CORE-BEG2007-PHASE_3"
                        # intake patterns
                        # At phase 3 ac is consumed
                        good_pattern = true
                        for exch in ["EX_ac_e"]
                            flx = _core_sol[RIDX[exch]]
                            if flx > 1e-2 
                                good_pattern = false
                                break
                            end
                        end
                        good_pattern || continue
                    end
                    
                    # --------------------------------------------
                    # distribution filters
                    # --------------------------------------------
                    # biomass distribution
                    ixd = RIDX["BIOMASS_Ecoli_core_w_GAM"]
                    _core_biom = _core_sol[ixd]
                    _core_biom in BIOMAS_RANGE || continue

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
        _dat = Dict("ensem" => ensem, "RIDX" => RIDX, "BIOMAS_RANGE" => BIOMAS_RANGE)
        fn = procdir(PROJ, 
            ["ensembles"], 
            SIMVER, ENSEM_SIZE,
            (;h = hash(ensem)), 
            ".jls"
        )
        sdat(_dat, fn; verbose = true)

        _ensem_summary(ensem, core_lep0)
        
        nothing

    end
end

