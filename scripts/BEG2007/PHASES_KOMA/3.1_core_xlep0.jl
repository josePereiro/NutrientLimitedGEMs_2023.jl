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
@tempcontext ["CORE_XLEP" => v"0.1.0"] let

    # -------------------------------------------
    # globals
    GLOB_DB = query(["ROOT", "GLOBALS"])
    LP_SOLVER = GLOB_DB["LP_SOLVER"]
    CORE_NET_ID = GLOB_DB["CORE_NET_ID"]
    GEM_NET_ID = GLOB_DB["GEM_NET_ID"]
    NTHREADS = GLOB_DB["NTHREADS"]

    # -------------------------------------------
    # nets
    # TODO: load 
    core_net0 = pull_net("ecoli_core")
    gem_net0 = pull_net("iJO1366")

    # -------------------------------------------
    # Add nutrient uptake reations
    core_net0 = resize(core_net0; 
        nmets = size(core_net0, 1) + 20,
        nrxns = size(core_net0, 2) + 20,
    )

    # -------------------------------------------
    # Galactose

    println("-"^40)
    println("Galactose")
    println("-"^40)

    # rxn[935]: EX_gal_e (D-Galactose exchange)
    # subsys: nothing
    # lb: 0.0, ub: 1000.0
    # (-1.0) gal_e ==> 
    colid = "EX_gal_e"
    @assert !hascolid(core_net0, colid)
    set_constraint!(core_net0, colid;
        S = Dict("gal_e" => -1),
        lb = 0.0, ub = 1000.0,
    )
    _merge_rxninfo!(core_net0, "EX_gal_e", gem_net0, "EX_gal_e")
    _merge_metinfo!(core_net0, "gal_e", gem_net0, "gal_e")
    println(colid, ": ", col_str(core_net0, colid))

    # rxn[1310]: GALtex (D-galactose transport via diffusion (extracellular to periplasm))
    # subsys: nothing
    # lb: -1000.0, ub: 1000.0
    # (-1.0) gal_e <==> (1.0) gal_p
    # Not apply

    # rxn[1308]: GALabcpp (D-galactose transport via ABC system (periplasm))
    # subsys: nothing
    # lb: 0.0, ub: 1000.0
    # (-1.0) atp_c + (-1.0) gal_p + (-1.0) h2o_c ==> (1.0) adp_c + (1.0) gal_c + (1.0) h_c + (1.0) pi_c
    colid = "GALabcpp"
    @assert !hascolid(core_net0, colid)
    @assert hasrowid(core_net0, "gal_e")
    @assert hasrowid(core_net0, "atp_c")
    @assert hasrowid(core_net0, "h2o_e")
    @assert hasrowid(core_net0, "adp_c")
    @assert !hasrowid(core_net0, "gal_c")
    @assert hasrowid(core_net0, "h_c")
    @assert hasrowid(core_net0, "pi_c")
    set_constraint!(core_net0, colid;
        S = Dict(
            "atp_c" => -1, "gal_e" => -1, "h2o_e" => -1, 
            "adp_c" => 1, "gal_c" => 1, "h_c" => 1, "pi_c" => 1, 
        ),
        lb = 0.0, ub = 1000.0,
    )
    _merge_rxninfo!(core_net0, "GALabcpp", gem_net0, "GALabcpp")
    _merge_metinfo!(core_net0, "gal_c", gem_net0, "gal_c")
    println(colid, ": ", col_str(core_net0, colid))

    # rxn[1299]: GALKr (Galactokinase)
    # subsys: nothing
    # lb: -1000.0, ub: 1000.0
    # (-1.0) atp_c + (-1.0) gal_c <==> (1.0) adp_c + (1.0) gal1p_c + (1.0) h_c
    colid = "GALKr"
    @assert !hascolid(core_net0, colid)
    @assert hasrowid(core_net0, "gal_c")
    @assert !hasrowid(core_net0, "gal1p_c")
    set_constraint!(core_net0, colid;
        S = Dict(
            "atp_c" => -1, "gal_c" => -1, 
            "adp_c" => 1, "gal1p_c" => 1, "h_c" => 1
        ),
        lb = -1000.0, ub = 1000.0,
    )
    _merge_rxninfo!(core_net0, "GALKr", gem_net0, "GALKr")
    _merge_metinfo!(core_net0, "gal1p_c", gem_net0, "gal1p_c")
    println(colid, ": ", col_str(core_net0, colid))

    # rxn[2522]: UGLT (UDPglucose--hexose-1-phosphate uridylyltransferase)
    # subsys: nothing
    # lb: -1000.0, ub: 1000.0
    # (-1.0) gal1p_c + (-1.0) udpg_c <==> (1.0) g1p_c + (1.0) udpgal_c
    colid = "UGLT"
    @assert !hascolid(core_net0, colid)
    @assert hasrowid(core_net0, "gal1p_c")
    @assert !hasrowid(core_net0, "udpg_c")
    @assert !hasrowid(core_net0, "g1p_c")
    @assert !hasrowid(core_net0, "udpgal_c")
    set_constraint!(core_net0, colid;
        S = Dict(
            "gal1p_c" => -1, "udpg_c" => -1, 
            "g1p_c" => 1, "udpgal_c" => 1
        ),
        lb = -1000.0, ub = 1000.0,
    )
    _merge_rxninfo!(core_net0, "UGLT", gem_net0, "UGLT")
    _merge_metinfo!(core_net0, "udpg_c", gem_net0, "udpg_c")
    _merge_metinfo!(core_net0, "g1p_c", gem_net0, "g1p_c")
    _merge_metinfo!(core_net0, "udpgal_c", gem_net0, "udpgal_c")
    println(colid, ": ", col_str(core_net0, colid))

    # rxn[2511]: UDPG4E (UDPglucose 4-epimerase)
    # subsys: nothing
    # lb: -1000.0, ub: 1000.0
    # (-1.0) udpg_c <==> (1.0) udpgal_c
    colid = "UDPG4E"
    @assert !hascolid(core_net0, colid)
    @assert hasrowid(core_net0, "udpg_c")
    @assert hasrowid(core_net0, "udpgal_c")
    @assert hasrowid(core_net0, "g1p_c")
    set_constraint!(core_net0, colid;
        S = Dict(
            "udpg_c" => -1, 
            "udpgal_c" => 1,
        ),
        lb = -1000.0, ub = 1000.0,
    )
    _merge_rxninfo!(core_net0, "UDPG4E", gem_net0, "UDPG4E")
    println(colid, ": ", col_str(core_net0, colid))

    # rxn[2082]: PGMT (Phosphoglucomutase)
    # subsys: nothing
    # lb: -1000.0, ub: 1000.0
    # (-1.0) g1p_c <==> (1.0) g6p_c
    colid = "PGMT"
    @assert !hascolid(core_net0, colid)
    @assert hasrowid(core_net0, "g6p_c")
    set_constraint!(core_net0, colid;
        S = Dict(
            "g1p_c" => -1, 
            "g6p_c" => 1,
        ),
        lb = -1000.0, ub = 1000.0,
    )
    _merge_rxninfo!(core_net0, "PGMT", gem_net0, "PGMT")
    println(colid, ": ", col_str(core_net0, colid))

    # -------------------------------------------
    # Maltose

    println()
    println("-"^40)
    println("Maltose")
    println("-"^40)

    # rxn[1000]: EX_malt_e (Maltose exchange)
    # subsys: nothing
    # lb: 0.0, ub: 1000.0
    # (-1.0) malt_e ==> 
    colid = "EX_malt_e"
    @assert !hascolid(core_net0, colid)
    @assert !hasrowid(core_net0, "malt_e")
    set_constraint!(core_net0, colid;
        S = Dict(
            "malt_e" => -1, 
        ),
        lb = -0.0, ub = 1000.0,
    )
    _merge_rxninfo!(core_net0, "EX_malt_e", gem_net0, "EX_malt_e")
    _merge_metinfo!(core_net0, "malt_e", gem_net0, "malt_e")
    println(colid, ": ", col_str(core_net0, colid))

    # rxn[1726]: MALTtexi (MaltoseMaltotriose transport via diffusion (extracellular to periplasm) irreversible)
    # subsys: nothing
    # lb: 0.0, ub: 1000.0
    # (-1.0) malt_e ==> (1.0) malt_p
    # Not apply

    # malEGFK
    # rxn[1724]: MALTabcpp (Maltose transport via ABC system (periplasm))
    # subsys: nothing
    # lb: 0.0, ub: 1000.0
    # (-1.0) atp_c + (-1.0) h2o_c + (-1.0) malt_p ==> (1.0) adp_c + (1.0) h_c + (1.0) malt_c + (1.0) pi_c
    colid = "MALTabcpp"
    @assert !hascolid(core_net0, colid)
    @assert hasrowid(core_net0, "malt_e")
    @assert !hasrowid(core_net0, "malt_c")
    set_constraint!(core_net0, colid;
        S = Dict(
            "atp_c" => -1, "h2o_c" => -1, "malt_e" => -1, 
            "adp_c" => 1, "h_c" => 1, "malt_c" => 1, "pi_c" => 1
        ),
        lb = -0.0, ub = 1000.0,
    )
    _merge_rxninfo!(core_net0, "MALTabcpp", gem_net0, "MALTabcpp")
    _merge_metinfo!(core_net0, "malt_c", gem_net0, "malt_c")
    println(colid, ": ", col_str(core_net0, colid))

    # NOTE
    # rxn[1725]: MALTptspp (Maltose transport via PEP:Pyr PTS (periplasm))
    # subsys: nothing
    # lb: 0.0, ub: 1000.0
    # (-1.0) malt_p + (-1.0) pep_c ==> (1.0) malt6p_c + (1.0) pyr_c
    # NOTE: malt6p_c is a dead end at iJO1366, so, using MALTabcpp

    # rxn[318]: AMALT1 (Amylomaltase (maltotriose))
    # subsys: nothing
    # lb: 0.0, ub: 1000.0
    # (-1.0) malt_c + (-1.0) malttr_c ==> (1.0) glc__D_c + (1.0) maltttr_c
    # I will ignore maltttr_c
    colid = "AMALT1"
    @assert !hascolid(core_net0, colid)
    @assert hasrowid(core_net0, "glc__D_e")
    @assert hasrowid(core_net0, "malt_c")
    set_constraint!(core_net0, colid;
        S = Dict(
            "malt_c" => -1,
            "glc__D_e" => 2,
        ),
        lb = -0.0, ub = 1000.0,
    )
    _merge_rxninfo!(core_net0, "AMALT1", gem_net0, "AMALT1")
    println(colid, ": ", col_str(core_net0, colid))

    # -------------------------------------------
    # Glycerol

    println()
    println("-"^40)
    println("Glycerol")
    println("-"^40)

    # rxn[958]: EX_glyc_e (Glycerol exchange)
    # subsys: nothing
    # lb: 0.0, ub: 1000.0
    # (-1.0) glyc_e ==> 
    colid = "EX_glyc_e"
    @assert !hascolid(core_net0, colid)
    @assert !hasrowid(core_net0, "glyc_e")
    set_constraint!(core_net0, colid;
        S = Dict(
            "glyc_e" => -1,
        ),
        lb = -0.0, ub = 1000.0,
    )
    _merge_rxninfo!(core_net0, "EX_glyc_e", gem_net0, "EX_glyc_e")
    _merge_metinfo!(core_net0, "glyc_e", gem_net0, "glyc_e")
    println(colid, ": ", col_str(core_net0, colid))

    # rxn[1406]: GLYCtex (Glycerol transport via diffusion (extracellular to periplasm))
    # subsys: nothing
    # lb: -1000.0, ub: 1000.0
    # (-1.0) glyc_e <==> (1.0) glyc_p
    # Not apply

    # rxn[1407]: GLYCtpp (Glycerol transport via channel (periplasm))
    # subsys: nothing
    # lb: -1000.0, ub: 1000.0
    # (-1.0) glyc_c <==> (1.0) glyc_p
    colid = "GLYCtpp"
    @assert !hascolid(core_net0, colid)
    @assert hasrowid(core_net0, "glyc_e")
    @assert !hasrowid(core_net0, "glyc_c")
    set_constraint!(core_net0, colid;
        S = Dict(
            "glyc_e" => -1,
            "glyc_c" => 1,
        ),
        lb = -1000.0, ub = 1000.0,
    )
    _merge_rxninfo!(core_net0, "GLYCtpp", gem_net0, "GLYCtpp")
    _merge_metinfo!(core_net0, "glyc_c", gem_net0, "glyc_c")
    println(colid, ": ", col_str(core_net0, colid))

    # rxn[1408]: GLYK (Glycerol kinase)
    # subsys: nothing
    # lb: 0.0, ub: 1000.0
    # (-1.0) atp_c + (-1.0) glyc_c ==> (1.0) adp_c + (1.0) glyc3p_c + (1.0) h_c
    colid = "GLYK"
    @assert !hascolid(core_net0, colid)
    @assert hasrowid(core_net0, "glyc_c")
    @assert !hasrowid(core_net0, "glyc3p_c")
    set_constraint!(core_net0, colid;
        S = Dict(
            "atp_c" => -1, "glyc_c" => -1,
            "adp_c" => 1, "glyc3p_c" => 1, "h_c" => 1
        ),
        lb = 0.0, ub = 1000.0,
    )
    _merge_rxninfo!(core_net0, "GLYK", gem_net0, "GLYK")
    _merge_metinfo!(core_net0, "glyc3p_c", gem_net0, "glyc3p_c")
    println(colid, ": ", col_str(core_net0, colid))

    # rxn[1267]: G3PD2 (Glycerol-3-phosphate dehydrogenase (NADP))
    # subsys: nothing
    # lb: -1000.0, ub: 1000.0
    # (-1.0) glyc3p_c + (-1.0) nadp_c <==> (1.0) dhap_c + (1.0) h_c + (1.0) nadph_c
    colid = "G3PD2"
    @assert !hascolid(core_net0, colid)
    @assert hasrowid(core_net0, "glyc3p_c")
    @assert hasrowid(core_net0, "dhap_c")
    @assert hasrowid(core_net0, "nadp_c")
    @assert hasrowid(core_net0, "nadph_c")
    set_constraint!(core_net0, colid;
        S = Dict(
            "glyc3p_c" => -1, "nadp_c" => -1,
            "dhap_c" => 1, "h_c" => 1, "nadph_c" => 1
        ),
        lb = -1000.0, ub = 1000.0,
    )
    _merge_rxninfo!(core_net0, "G3PD2", gem_net0, "G3PD2")
    println(colid, ": ", col_str(core_net0, colid))

    # -------------------------------------------
    # open complex medium
    core_net0 = emptyless_model(core_net0)

    # Bounds from Beg et al, 2007, fig 3
    # We just need the maximum rates. 
    # Original in (mmol/ min g)
    # 1 [mmol/ min g] * 60 = 1 [mmol/ h g]
    # adjustment growth factor (I aim to adjust the maximum growth to match glucose only regime ~0.8 1/h (Fig2, a))
    # At least I keep the relation between maximal nutrient intakes
    @stage! gf = 0.18 
    if SIMVER == "ECOLI-CORE-BEG2007-PHASE_0"
        lb!(core_net0, "EX_glc__D_e", -0.9 * 60 * gf)
        lb!(core_net0, "EX_lac__D_e", -1.0 * 60 * gf)
        lb!(core_net0, "EX_malt_e", -0.1 * 60 * gf)
        lb!(core_net0, "EX_gal_e", -0.2 * 60 * gf)
        lb!(core_net0, "EX_glyc_e", -0.6 * 60 * gf)
        lb!(core_net0, "EX_ac_e", 0.0)
        ub!(core_net0, "EX_ac_e", 2.0 * 60 * gf)
    elseif SIMVER == "ECOLI-CORE-BEG2007-PHASE_1"
        lb!(core_net0, "EX_glc__D_e", -0.9 * 60 * gf)
        lb!(core_net0, "EX_lac__D_e", 0.0)
        lb!(core_net0, "EX_malt_e", 0.0)
        lb!(core_net0, "EX_gal_e", 0.0)
        lb!(core_net0, "EX_glyc_e", 0.0)
        lb!(core_net0, "EX_ac_e", 0.0)
        ub!(core_net0, "EX_ac_e", 2.0 * 60 * gf)
    elseif SIMVER == "ECOLI-CORE-BEG2007-PHASE_2"
        lb!(core_net0, "EX_glc__D_e", 0.0)
        lb!(core_net0, "EX_lac__D_e", -1.0 * 60 * gf)
        lb!(core_net0, "EX_malt_e", -0.1 * 60 * gf)
        lb!(core_net0, "EX_gal_e", -0.2 * 60 * gf)
        lb!(core_net0, "EX_glyc_e", -0.6 * 60 * gf)
        lb!(core_net0, "EX_ac_e", -1.5 * 60 * gf) # (?)
        ub!(core_net0, "EX_ac_e", 2.0 * 60 * gf)
    elseif SIMVER == "ECOLI-CORE-BEG2007-PHASE_3"
        lb!(core_net0, "EX_glc__D_e", 0.0)
        lb!(core_net0, "EX_lac__D_e", 0.0)
        lb!(core_net0, "EX_malt_e", 0.0)
        lb!(core_net0, "EX_gal_e", 0.0)
        lb!(core_net0, "EX_glyc_e", -0.6 * 60 * gf)
        lb!(core_net0, "EX_ac_e", -1.5 * 60 * gf)
        ub!(core_net0, "EX_ac_e", 2.0 * 60 * gf)
    end

    # -------------------------------------------
    # Post processing
    core_lep0 = lepmodel(core_net0)
    
    # blep
    cid = (:BOX, CORE_NET_ID, hash(core_lep0))
    _, core_blep0ref = withcachedat(PROJ, :get!, cid) do 
        core_blep0 = box(core_lep0, LP_SOLVER; nths = NTHREADS, verbose = true)
        return CacheRef(core_blep0)
    end

    # elep
    cid = (:ELEP, CORE_NET_ID, hash(core_blep0ref))
    _, core_elep0ref = withcachedat(PROJ, :get!, cid) do 
        core_blep0 = core_blep0ref[]
        core_elep0 = EchelonLEPModel(core_blep0; verbose = true)
        return CacheRef(core_elep0)
    end

    # stage nets
    core_lep0ref = CacheRef(core_lep0)
    core_net0ref = CacheRef(core_net0)
    @stage! "core_net0" => core_net0ref
    @stage! "core_lep0" => core_lep0ref
    @stage! "core_blep0" => core_blep0ref
    @stage! "core_elep0" => core_elep0ref

    # -------------------------------------------
    # FBA TEST

    println()
    println("="^40)
    println("LEP FBA TEST")
    println()

    let 
        _core_net0 = deepcopy(core_net0)
        biom_id = extras(_core_net0, "BIOM")
        linear_weights!(_core_net0, biom_id, 1.0)
        exchs = ["EX_glc__D_e", "EX_lac__D_e", "EX_malt_e", "EX_gal_e", "EX_glyc_e", "EX_ac_e"]
        lb!(_core_net0, exchs, 0.0)

        for nut_id in exchs

            println("-"^40)
            println(nut_id)

            lb!(_core_net0, nut_id, -10.0)
            opm = fba(_core_net0, Clp.Optimizer)
            println(biom_id, ": ", solution(opm, biom_id))
            # println(nut_id, ": ", solution(opm, nut_id))
            @assert solution(opm, biom_id) > 1e-2
            lb!(_core_net0, nut_id, 0.0)
        end    
    end

    # Test FBA
    println()
    println("="^40)
    println("ELEP FBA TEST")
    println()
    
    let 
        _core_lep0 = core_lep0ref[]
        _core_elep0 = core_elep0ref[]
        obj_id = extras(_core_lep0, "BIOM")
        _opm = FBAOpModel(_core_lep0, LP_SOLVER)
        set_linear_obj!(_opm, obj_id, MAX_SENSE)
        optimize!(_opm)
        @show solution(_opm, obj_id)
        @assert solution(_opm, obj_id) > 0.5
        
        @show length(_core_elep0.idxi)
        @show size(_core_lep0, 2)
        @show length(_core_elep0.idxi) / size(_core_lep0, 2)
    end

    nothing
end

## ------------------------------------------------------------
# save
_save_contextdb(SIMVER)