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
# ------------------------------------------------------------
include("1_setup_sim.jl")
include("1.1_utils.jl")

## ------------------------------------------------------------
let
    # koma files
    downreg_factor = 0.3 # TOSYNC
    objfiles = readdir(procdir(PROJ, [SIMVER]); join = true)
    @threads for fn in sort(objfiles)
        contains(basename(fn), "obj_reg") || continue

        # deserialize
        _, obj_reg = ldat(fn)
        isempty(obj_reg) && continue
        
        info_frec = 1000
        do_save = false
        for (obji, obj) in enumerate(obj_reg)

            # filter
            haskey(obj, "core_fva.ver") || continue
            do_save = true

            # info
            show_flag = obji == 1 || obji == lastindex(obj_reg) || iszero(rem(obji, info_frec)) 
            show_flag && println("[", getpid(), ".", threadid(), "] ", 
                "obji ", obji, "\\", length(obj_reg), " ",
                basename(fn)
            )

            # format
            # DEPRECATE
            _DEP = get!(obj, "_DEPRECATED", Dict())
            _DEP["v1"] = Dict()
            _DEP["v1"]["core_fva.desc"] = "fva of strip.koset[1:end-1]"
            _DEP["v1"]["core_fva.fvalb"] = obj["core_fva.fvalb"]
            _DEP["v1"]["core_fva.fvaub"] = obj["core_fva.fvaub"]
            _DEP["v1"]["core_fva.kos"] = obj["core_fva.kos"]
            _DEP["v1"]["core_fva.ver"] = obj["core_fva.ver"]

            # DELETE
            delete!(obj, "core_fva.fvalb")
            delete!(obj, "core_fva.fvaub")
            delete!(obj, "core_fva.kos")
            delete!(obj, "core_fva.ver")

        end # for obj
        
        # return obj_reg
        do_save && sdat(obj_reg, fn)

    end # for fn 

end