## ------------------------------------------------------------
@time begin
    using MetXGEMs
    using MetXBase
    using MetXOptim
    using MetXNetHub
    using NutrientLimitedGEMs
end

# ------------------------------------------------------------
include("1_setup_sim.jl")
include("1.1_utils.jl")

## ------------------------------------------------------------
# Prepare network
@tempcontext ["CORE_XLEP" => v"0.1.0"] let

    # globals
    glob_db = query(["ROOT", "GLOBALS"])
    LP_SOLVER = glob_db["LP_SOLVER"]
    CORE_NET_ID = glob_db["CORE_NET_ID"]
    NTHREADS = glob_db["NTHREADS"]

    # net
    net0 = pull_net(CORE_NET_ID)
    
    #  open complex medium
    # Bounds from Beg et al, 2007, fig 3
    # We just need the maximum rates. 
    # Original in (mmol/ min g)
    # 1 [mmol/ min g] * 60 = 1 [mmol/ h g]
    # adjustment growth factor (I aim to adjust the maximum growth to match glucose only regime ~0.8 1/h (Fig2, a))
    # At least I keep the relation between maximal nutrient intakes
    gf = 0.18 
    lb!(net0, "EX_glc__D_e", -0.9 * 60 * gf)
    lb!(net0, "EX_lac__D_e", -1.0 * 60 * gf)
    # # lb!(net0, "EX_malt_e", -0.1 * 60 * gf)
    # # lb!(net0, "EX_gal_e", -0.2 * 60 * gf)
    # # lb!(net0, "EX_glyc_e", -0.6 * 60 * gf)
    lb!(net0, "EX_ac_e", -1.5 * 60 * gf)
    ub!(net0, "EX_ac_e", 2.0 * 60 * gf)

    # opm = fba(net0, Clp.Optimizer)
    # @show solution(opm, extras(net0, "BIOM"))
    # @show solution(opm, extras(net0, "EX_GLC"))

    # lep
    lep0 = lepmodel(net0)
    
    # blep
    cid = (:BOX, hash(lep0))
    _, blep0ref = withcachedat(PROJ, :get!, cid) do 
        blep0 = box(lep0, LP_SOLVER; nths = NTHREADS, verbose = true)
        return CacheRef(blep0)
    end

    # elep
    cid = (:ELEP, CORE_NET_ID, hash(blep0ref))
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