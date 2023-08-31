## ------------------------------------------------------------
@time begin
    using Dates
    using Random
    using MetXGEMs
    using MetXBase
    using MetXOptim
    using MetXNetHub
    using Base.Threads
    using ProgressMeter
    using SimpleLockFiles
    using NutrientLimitedGEMs
end

# ------------------------------------------------------------
# TODO: refactor at koma.jl
let
    komafiles = readdir(procdir(PROJ, [SIMVER]); join = true)
    @threads for koma_fn in komafiles
        contains(basename(koma_fn), "koma_reg") || continue
        obj_fn = replace(koma_fn, "koma_reg" => "obj_reg")
        isfile(obj_fn) && continue
        @show obj_fn
        
        # deserialize
        _, koma_reg = ldat(koma_fn)
        obj_reg = Dict{String, Any}[]

        # reformat
        # TODO: fix this in koma script
        for (koset, st) in koma_reg
            obj = Dict(
                "koset" => koset,
                "koma_status" => st::Symbol,
            )
            push!(obj_reg, obj)
        end
        empty!(koma_reg)
        
        sdat(obj_reg, obj_fn)
        rm(koma_fn; force = true)
    end
end