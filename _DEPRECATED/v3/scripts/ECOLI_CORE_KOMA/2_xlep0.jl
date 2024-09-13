## ------------------------------------------------------------
@time begin
    using MetXGEMs
    using MetXBase
    using MetXOptim
    using MetXNetHub
    using NutrientLimitedGEMs_2023
end

# ------------------------------------------------------------
include("1_setup_sim.jl")
include("1.1_utils.jl")

## ------------------------------------------------------------
# Prepare network
@tempcontext ["XLEP" => v"0.1.0"] let

    # globals
    glob_db = query(["ROOT", "GLOBALS"])
    LP_SOLVER = glob_db["LP_SOLVER"]
    NET_ID = glob_db["NET_ID"]
    NTHREADS = glob_db["NTHREADS"]

    # net
    net0 = pull_net(NET_ID)
    
    #  open complex medium
    lb!(net0, "EX_glc__D_e", -10.0)
    lb!(net0, "EX_lac__D_e", -10.0)
    # lb!(net0, "EX_malt_e", -10.0)
    # lb!(net0, "EX_gal_e", -10.0)
    # lb!(net0, "EX_glyc_e", -10.0)
    lb!(net0, "EX_ac_e", -10.0)

    # lep
    lep0 = lepmodel(net0)
    
    # blep
    cid = (:BOX, hash(lep0))
    _, blep0ref = withcachedat(PROJ, :get!, cid) do 
        blep0 = box(lep0, LP_SOLVER; nths = NTHREADS, verbose = true)
        return CacheRef(blep0)
    end

    # elep
    cid = (:ELEP, NET_ID, hash(blep0ref))
    _, elep0ref = withcachedat(PROJ, :get!, cid) do 
        blep0 = blep0ref[]
        elep0 = EchelonLEPModel(blep0; verbose = true)
        return CacheRef(elep0)
    end

    @stage! "net0" => CacheRef(net0)
    @stage! "lep0" => CacheRef(lep0)
    @stage! "blep0" => blep0ref
    @stage! "elep0" => elep0ref

    # Test FBA
    lep = lepmodel(elep0ref[])
    obj_id = extras(lep, "BIOM")
    @show obj_id
    opm = FBAOpModel(lep, LP_SOLVER)
    set_linear_obj!(opm, obj_id, MAX_SENSE)
    optimize!(opm)
    @show solution(opm, obj_id)
    @assert solution(opm, obj_id) > 0.5
    
    @show size(lep0, 2)
    @show length(elep0ref[].idxi)
    @show length(elep0ref[].idxi) / size(lep0, 2)

    nothing
end

## ------------------------------------------------------------
# save
_save_contextdb(SIMVER)