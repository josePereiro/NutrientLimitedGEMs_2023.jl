## ------------------------------------------------------------
@time begin
    using Clp
    using Gurobi
    using MetXGEMs
    using MetXBase
    using MetXOptim
    using MetXNetHub
    using NutrientLimitedGEMs
end

# ------------------------------------------------------------
include("1_setup_sim.jl")

## ------------------------------------------------------------
# Prepare network
# TODO: Do this and store it at MetXNetHub (document the source version)
@tempcontext ["ELEP" => v"0.1.0"] let
    # @stage! netid = "ecoli_core"
    @stage!  netid = "iJR904"
    # @stage! netid = "iCHO2291"
    # @stage! netid = "SysBioChalmers_Human_GEM"

    net0 = pull_net(netid)
    lep0 = lepmodel(net0)

    
    # @time global elep0 = EchelonLEPModel(lep0; verbose = true)
    cid = (:ELEP, netid, hash(lep0))
    _, elep0ref = withcachedat(PROJ, :get!, cid) do 
        
        # elep0 = EchelonLEPModel(lep0; verbose = true)
        # ToDel
        fn = joinpath(@__DIR__, string("elep0-", netid, ".jls"))
        elep0 = deserialize(fn)
        
        return CacheRef(elep0)
    end

    @stage! "net0" => CacheRef(net0)
    @stage! "lep0" => CacheRef(lep0)
    @stage! "elep0" => elep0ref
    
    @show size(lep0, 2)
    @show length(elep0ref[].idxi)
    @show length(elep0ref[].idxi) / size(lep0, 2)

    nothing
end

## ------------------------------------------------------------