## ------------------------------------------------------------
@time begin
    using MetXGEMs
    using MetXBase
    using MetXOptim
    using MetXNetHub
    using NutrientLimitedGEMs_2023
end

# ------------------------------------------------------------
include("1_setup.jl")

## ------------------------------------------------------------
# Prepare network
@tempcontext ["GEM_XLEP" => v"0.1.0"] let

    # globals
    GLOB_DB = query(["ROOT", "GLOBALS"])
    LP_SOLVER = GLOB_DB["LP_SOLVER"]
    GEM_NET_ID = GLOB_DB["GEM_NET_ID"]
    NTHREADS = GLOB_DB["NTHREADS"]

    # net
    net0 = pull_net(GEM_NET_ID)
    
    #  open complex medium
    lb!(net0, "R_EX_glc__D_e", -10.0)
    lb!(net0, "R_EX_lac__L_e", -10.0)
    lb!(net0, "R_EX_malt_e", -10.0)
    lb!(net0, "R_EX_gal_e", -10.0)
    lb!(net0, "R_EX_glyc_e", -10.0)
    lb!(net0, "R_EX_ac_e", -10.0)

    # lep
    lep0 = lepmodel(net0)
    
    # blep
    cid = (:BOX, hash(lep0))
    _, blep0ref = withcachedat(PROJ, :get!, cid) do 
        blep0 = box(lep0, LP_SOLVER; nths = NTHREADS, verbose = true)
        return CacheRef(blep0)
    end

    # elep
    cid = (:ELEP, GEM_NET_ID, hash(blep0ref))
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