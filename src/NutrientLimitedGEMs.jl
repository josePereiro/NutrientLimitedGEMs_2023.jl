module NutrientLimitedGEMs

    using ProjAssistant
    using Gurobi
    using MetX
    using Plots

    include("utils.jl")
    
    @gen_top_proj()

    function __init__()
        @create_proj_dirs()
    end

end
