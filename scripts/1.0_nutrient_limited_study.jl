## ------------------------------------------------------------------
# Launch workers
using Distributed
while length(workers()) < 3
    addprocs(1; 
        exeflags = ["--startup-file=no", "--project=$(Base.active_project())"]
    )
end

## ------------------------------------------------------------------
@everywhere function _do_import()
    @eval Main begin
        using COBREXA
        using NutrientLimitedGEMs
        const NL = NutrientLimitedGEMs
        using GLPK
        using ProjAssistant
    end
end

## ------------------------------------------------------------------
# Precompile first
@time _do_import()

## ------------------------------------------------------------------
# import everywhere
@everywhere begin
    println("Active project: ", Base.active_project())
    @time _do_import()
end

@info("Working at $(workers())")

## ------------------------------------------------------------------
function _exchs_compute_effect_mat(model, nsample; 
        abs_max = 100.0, kwargs...
    )
    
    # compute dep mat
    for exchid in find_uptake_exchs(model)
        println("\n", "-"^60)

        # samples
        lb, ub = bounds(model, exchid)
        lb, ub = max(lb, -abs_max), min(abs(lb), abs_max)
        samples = lb .+ (range(0.0, 1.0; length = nsample) .* (ub - lb))
        
        # compute
        @time vals, mat = compute_effect_mat(
            model, exchid, samples; 
            kwargs...
        )
    end

    @info("DONE")
end

## ------------------------------------------------------------------
# HartGEMs
# medium dep
let

    return

    # load model
    model_params = (;modelid="stst_scaled", stst="E",  tissue="GBM")
    # model_params = (;modelid="fva_scaled",  tissue="GBM")
    # model_params = (;modelid="fva_base",  tissue="GBM")
    @time model = load_HartGEM_model(; model_params...)
    nsample = 50

    # save globals
    sdat(NL, 
        (;model_params, nsample), 
        "effect_mat-globals-HartGEMs.jls"
    )
    
    _exchs_compute_effect_mat(model, nsample; 
        abs_max = 100.0, 
        workers = workers(),
        reltol = 0.01, 
    )

end


## ------------------------------------------------------------------
# EColi
# medium dep
let
    return

    # load model
    model_params = (;modelid="ecolicore")
    @time model = load_ecoli()
    nsample = 500

    # save globals
    sdat(NL, 
        (;model_params, nsample), 
        "effect_mat-globals-EColi.jls"
    )
    
    _exchs_compute_effect_mat(model, nsample; 
        clearcache = false, 
        abs_max = 100.0
    )

end

## ------------------------------------------------------------------
# Toy models
# medium dep
let

    # load model
    model_params = (;modelid="toy_model3")
    @time model = NL.load_toy_model3()
    nsample = 500

    # save globals
    sdat(NL, 
        (;model_params, nsample), 
        "effect_mat-globals-ToyModel.jls"
    )
    
    _exchs_compute_effect_mat(model, nsample; 
        clearcache = true, 
        abs_max = 100.0
    )

end