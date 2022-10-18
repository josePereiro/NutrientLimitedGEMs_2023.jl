module NutrientLimitedGEMs

    using MetNets, MetLP, MetEP
    using ProjAssistant
    using Distributed
    using Serialization
    import GLPK
    using Plots
    import OrderedCollections
    import OrderedCollections: OrderedDict
    using Base.Threads
    using Distributions
    
    @gen_top_proj()

    include("Utils/utils.jl")
    include("Utils/linear_fit.jl")
    include("Utils/lineal_dep.jl")
    include("Utils/connectome.jl")
    include("Utils/fixxed_mat.jl")
    include("Utils/entropy.jl")
    include("Utils/delta_mat.jl")

    include("ToyModel/toy_model.jl")
    include("ToyModel/toy_model2.jl")
    include("ToyModel/toy_model3.jl")
    include("ToyModel/toy_model4.jl")
    include("ToyModel/toy_model5.jl")
    
    include("iJR904/load_iJR904_model.jl")
    
    include("EColiCore/EColiCore.jl")
    
    include("HartGEM/HartGEM.jl")

    # export everything that isn't prefixed with _ (inspired by JuMP.jl, thanks!)
    for sym in names(@__MODULE__, all = true)
        sym in [Symbol(@__MODULE__), :eval, :include] && continue
        startswith(string(sym), ['_', '#']) && continue
        @eval export $sym
    end

    function __init__()
        @create_proj_dirs()
    end

end
