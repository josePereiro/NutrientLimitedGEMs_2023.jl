## ------------------------------------------------------------
@time begin
    using MetXGEMs
    using MetXBase
    using MetXOptim
    using MetXNetHub
    using NutrientLimitedGEMs
end

# ------------------------------------------------------------
include("1_setup.jl")
include("2_utils.jl")

## ------------------------------------------------------------
# Prepare network
@tempcontext ["GEM_XLEP" => v"0.1.0"] let

    # -------------------------------------------
    # globals
    GLOB_DB = query(["ROOT", "GLOBALS"])
    LP_SOLVER = GLOB_DB["LP_SOLVER"]
    GEM_NET_ID = GLOB_DB["GEM_NET_ID"]
    NTHREADS = GLOB_DB["NTHREADS"]

    # -------------------------------------------
    # core_xlepdb
    core_xlepdb = query("ROOT", "CORE_XLEP")

    # -------------------------------------------
    # nets
    core_net0 = core_xlepdb["core_net0"][]
    gem_net0 = pull_net("iJO1366")

    # -------------------------------------------
    # open complex medium

    # Bounds from Beg et al, 2007, fig 3
    # see core xlep script
    # @stage! gf = 0.18 
    exch = "EX_glc__D_e"
    lb!(gem_net0, exch, lb(core_net0, exch))
    exch = "EX_lac__D_e"
    lb!(gem_net0, exch, lb(core_net0, exch))
    exch = "EX_malt_e"
    lb!(gem_net0, exch, lb(core_net0, exch))
    exch = "EX_gal_e"
    lb!(gem_net0, exch, lb(core_net0, exch))
    exch = "EX_glyc_e"
    lb!(gem_net0, exch, lb(core_net0, exch))
    exch = "EX_ac_e"
    lb!(gem_net0, exch, lb(core_net0, exch))
    exch = "EX_ac_e"
    ub!(gem_net0, exch, ub(core_net0, exch))

    # -------------------------------------------
    # Post processing
    gem_lep0 = lepmodel(gem_net0)
    
    # blep
    cid = (:BOX, GEM_NET_ID, hash(gem_lep0))
    _, gem_blep0ref = withcachedat(PROJ, :get!, cid) do 
        gem_blep0 = box(gem_lep0, LP_SOLVER; nths = NTHREADS, verbose = true)
        return CacheRef(gem_blep0)
    end

    # elep
    cid = (:ELEP, GEM_NET_ID, hash(gem_blep0ref))
    _, gem_elep0ref = withcachedat(PROJ, :get!, cid) do 
        gem_blep0 = gem_blep0ref[]
        gem_elep0 = EchelonLEPModel(gem_blep0; verbose = true)
        return CacheRef(gem_elep0)
    end

    # stage nets
    gem_lep0ref = CacheRef(gem_lep0)
    gem_net0ref = CacheRef(gem_net0)
    @stage! "gem_net0" => gem_net0ref
    @stage! "gem_lep0" => gem_lep0ref
    @stage! "gem_blep0" => gem_blep0ref
    @stage! "gem_elep0" => gem_elep0ref

    # -------------------------------------------
    # FBA TEST

    println()
    println("="^40)
    println("LEP FBA TEST")
    println()

    let 
        _gem_net0 = deepcopy(gem_net0)
        biom_id = extras(_gem_net0, "BIOM")
        linear_weights!(_gem_net0, biom_id, 1.0)
        exchs = ["EX_glc__D_e", "EX_lac__D_e", "EX_malt_e", "EX_gal_e", "EX_glyc_e", "EX_ac_e"]
        lb!(_gem_net0, exchs, 0.0)

        for nut_id in exchs

            println("-"^40)
            println(nut_id)

            lb!(_gem_net0, nut_id, -10.0)
            opm = fba(_gem_net0, Clp.Optimizer)
            println(biom_id, ": ", solution(opm, biom_id))
            # println(nut_id, ": ", solution(opm, nut_id))
            @assert solution(opm, biom_id) > 1e-2
            lb!(_gem_net0, nut_id, 0.0)
        end    
    end

    # Test FBA
    println()
    println("="^40)
    println("ELEP FBA TEST")
    println()
    
    let 
        _gem_lep0 = gem_lep0ref[]
        _gem_elep0 = gem_elep0ref[]
        obj_id = extras(_gem_lep0, "BIOM")
        _opm = FBAOpModel(_gem_lep0, LP_SOLVER)
        set_linear_obj!(_opm, obj_id, MAX_SENSE)
        optimize!(_opm)
        @show solution(_opm, obj_id)
        @assert solution(_opm, obj_id) > 0.5
        
        @show length(_gem_elep0.idxi)
        @show size(_gem_lep0, 2)
        @show length(_gem_elep0.idxi) / size(_gem_lep0, 2)
    end

    nothing
end

## ------------------------------------------------------------
# save
_save_contextdb(SIMVER)