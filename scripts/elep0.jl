## ------------------------------------------------------------
@time begin
    using MetXBase
    using MetXOptim
    using MetXNetHub
    using MetXGEMs
    using Base.Threads
    using Serialization
end

## ------------------------------------------------------------
# Prepare network
let
    # netid = "ecoli_core"
    # netid = "iJR904"
    netid = "SysBioChalmers_Human_GEM"
    @show netid
    lep0 = lepmodel(pull_net(netid))
    @time elep0 = EchelonLEPModel(lep0; verbose = true)
    @show size(lep0, 2)

    fn = joinpath(@__DIR__, string("elep0-",netid, ".jls"))
    serialize(fn, elep0)

    for ext in [".xml", ".mat"]
        fn = joinpath(@__DIR__, string("elep0-",netid, ext))
        @show fn
        save_net(elep0, fn)
    end

    nothing
end
