## ------------------------------------------------------------
@time begin
    using Random
    using MetXGEMs
    using MetXNetHub
    using MetXBase
    using MetXOptim
    using Base.Threads
    using NutrientLimitedGEMs
end

# TODO: Add Graphs.jl kind of functionality for getting basic stuff
# --------------------------------------------------------------
include("1_setup_sim.jl")
include("1.1_utils.jl")

## -------------------------------------------------------------
# I'll separate the batch files into single obje files
let
    bashsdir = procdir(PROJ, [SIMVER])
    objsdir = procdir(PROJ, [SIMVER, "objs"])
    mkdir(objsdir)
    bashfiles = readdir(bashsdir; join = true)
    # @threads 
    for bashfn in sort(bashfiles)
        contains(basename(bashfn), "obj_reg") || continue

    end
end

## -------------------------------------------------------------
let
    # koma files
    downreg_factor = 0.3 # TOSYNC
    objfiles = readdir(procdir(PROJ, [SIMVER]); join = true)
    # @threads 
    for fn in sort(objfiles)
        contains(basename(fn), "obj_reg") || continue

        # deserialize
        obj_reg = try_ldat(fn)
        isempty(obj_reg) && continue
        
        info_frec = 1000
        do_save = false
        for (obji, obj) in enumerate(obj_reg)

            # filter
            haskey(obj, "core_fva.ver") || continue
            haskey(obj, "core_biomass.ver") || continue
            haskey(obj, "core_nut_sp.ver") || continue
            haskey(obj, "core_feasets.ver") || continue
            haskey(obj, "core_feasets") || continue

            return obj

            # info
            show_flag = obji == 1 || obji == lastindex(obj_reg) || iszero(rem(obji, info_frec)) 
            show_flag && println("[", getpid(), ".", threadid(), "] ", 
                "obji ", obji, "\\", length(obj_reg), " ",
                basename(fn)
            )

            # format
            # DEPRECATE

            # DELETE
            for (li, feaobj) in obj["feasibles"]
                for key in feaobj
                    endswith(key, "vars_ms") || endswith(key, "vars_errs") || continue
                    feaobj[key] isa Vector{Float16} && continue
                    feaobj[key] = Float16.(feaobj[key])
                    do_save = true
                end
            end

        end # for obj
        
        # return obj_reg
        do_save && sdat(obj_reg, fn)

    end # for fn 

end