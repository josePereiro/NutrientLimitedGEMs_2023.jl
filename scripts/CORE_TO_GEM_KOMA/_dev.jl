## ------------------------------------------------------------
@time begin
    using MetXGEMs
    using MetXBase
    using MetXOptim
    using MetXNetHub
    using NutrientLimitedGEMs
    using Clp
    # using COBREXA
end

## ------------------------------------------------------------
let

    global cnet0 = pull_net("ecoli_core")
    global lep0 = resize(lepmodel(cnet0); 
        nrows = size(cnet0, 1) + 20,
        ncols = size(cnet0, 2) + 20,
    )
    # global gnet0 = pull_net("iJO1366")

    # -------------------------------------------
    # Galactose
    # -------------------------------------------
    println("-"^40)
    println("Galactose")
    println("-"^40)

    # rxn[935]: R_EX_gal_e (D-Galactose exchange)
    # subsys: nothing
    # lb: 0.0, ub: 1000.0
    # (-1.0) M_gal_e ==> 
    colid = "EX_gal_e"
    @assert !hascolid(lep0, colid)
    set_constraint!(lep0, colid;
        S = Dict("gal_e" => -1),
        lb = 0.0, ub = 1000.0,
    )
    println(colid, ": ", col_str(lep0, colid))

    # rxn[1310]: R_GALtex (D-galactose transport via diffusion (extracellular to periplasm))
    # subsys: nothing
    # lb: -1000.0, ub: 1000.0
    # (-1.0) M_gal_e <==> (1.0) M_gal_p
    # Not apply

    # rxn[1308]: R_GALabcpp (D-galactose transport via ABC system (periplasm))
    # subsys: nothing
    # lb: 0.0, ub: 1000.0
    # (-1.0) M_atp_c + (-1.0) M_gal_p + (-1.0) M_h2o_c ==> (1.0) M_adp_c + (1.0) M_gal_c + (1.0) M_h_c + (1.0) M_pi_c
    colid = "GALabcpp"
    @assert !hascolid(lep0, colid)
    @assert hasrowid(lep0, "gal_e")
    @assert hasrowid(lep0, "atp_c")
    @assert hasrowid(lep0, "h2o_e")
    @assert hasrowid(lep0, "adp_c")
    @assert !hasrowid(lep0, "gal_c")
    @assert hasrowid(lep0, "h_c")
    @assert hasrowid(lep0, "pi_c")
    set_constraint!(lep0, colid;
        S = Dict(
            "atp_c" => -1, "gal_e" => -1, "h2o_e" => -1, 
            "adp_c" => 1, "gal_c" => 1, "h_c" => 1, "pi_c" => 1, 
        ),
        lb = 0.0, ub = 1000.0,
    )
    println(colid, ": ", col_str(lep0, colid))

    # rxn[1299]: R_GALKr (Galactokinase)
    # subsys: nothing
    # lb: -1000.0, ub: 1000.0
    # (-1.0) M_atp_c + (-1.0) M_gal_c <==> (1.0) M_adp_c + (1.0) M_gal1p_c + (1.0) M_h_c
    colid = "GALKr"
    @assert !hascolid(lep0, colid)
    @assert hasrowid(lep0, "gal_c")
    @assert !hasrowid(lep0, "gal1p_c")
    set_constraint!(lep0, colid;
        S = Dict(
            "atp_c" => -1, "gal_c" => -1, 
            "adp_c" => 1, "gal1p_c" => 1, "h_c" => 1
        ),
        lb = -1000.0, ub = 1000.0,
    )
    println(colid, ": ", col_str(lep0, colid))

    # rxn[2522]: R_UGLT (UDPglucose--hexose-1-phosphate uridylyltransferase)
    # subsys: nothing
    # lb: -1000.0, ub: 1000.0
    # (-1.0) M_gal1p_c + (-1.0) M_udpg_c <==> (1.0) M_g1p_c + (1.0) M_udpgal_c
    colid = "UGLT"
    @assert !hascolid(lep0, colid)
    @assert hasrowid(lep0, "gal1p_c")
    @assert !hasrowid(lep0, "udpg_c")
    @assert !hasrowid(lep0, "g1p_c")
    @assert !hasrowid(lep0, "udpgal_c")
    set_constraint!(lep0, colid;
        S = Dict(
            "gal1p_c" => -1, "udpg_c" => -1, 
            "g1p_c" => 1, "udpgal_c" => 1
        ),
        lb = -1000.0, ub = 1000.0,
    )
    println(colid, ": ", col_str(lep0, colid))

    # rxn[2511]: R_UDPG4E (UDPglucose 4-epimerase)
    # subsys: nothing
    # lb: -1000.0, ub: 1000.0
    # (-1.0) M_udpg_c <==> (1.0) M_udpgal_c
    colid = "UDPG4E"
    @assert !hascolid(lep0, colid)
    @assert hasrowid(lep0, "udpg_c")
    @assert hasrowid(lep0, "udpgal_c")
    @assert hasrowid(lep0, "g1p_c")
    set_constraint!(lep0, colid;
        S = Dict(
            "udpg_c" => -1, 
            "udpgal_c" => 1,
        ),
        lb = -1000.0, ub = 1000.0,
    )
    println(colid, ": ", col_str(lep0, colid))

    # rxn[2082]: R_PGMT (Phosphoglucomutase)
    # subsys: nothing
    # lb: -1000.0, ub: 1000.0
    # (-1.0) M_g1p_c <==> (1.0) M_g6p_c
    colid = "PGMT"
    @assert !hascolid(lep0, colid)
    @assert hasrowid(lep0, "g6p_c")
    set_constraint!(lep0, colid;
        S = Dict(
            "g1p_c" => -1, 
            "g6p_c" => 1,
        ),
        lb = -1000.0, ub = 1000.0,
    )
    println(colid, ": ", col_str(lep0, colid))


    # -------------------------------------------
    # Maltose
    # -------------------------------------------
    println()
    println("-"^40)
    println("Maltose")
    println("-"^40)

    # rxn[1000]: R_EX_malt_e (Maltose exchange)
    # subsys: nothing
    # lb: 0.0, ub: 1000.0
    # (-1.0) M_malt_e ==> 
    colid = "EX_malt_e"
    @assert !hascolid(lep0, colid)
    @assert !hasrowid(lep0, "malt_e")
    set_constraint!(lep0, colid;
        S = Dict(
            "malt_e" => -1, 
        ),
        lb = -0.0, ub = 1000.0,
    )
    println(colid, ": ", col_str(lep0, colid))

    # rxn[1726]: R_MALTtexi (MaltoseMaltotriose transport via diffusion (extracellular to periplasm) irreversible)
    # subsys: nothing
    # lb: 0.0, ub: 1000.0
    # (-1.0) M_malt_e ==> (1.0) M_malt_p
    # Not apply

    # malEGFK
    # rxn[1724]: R_MALTabcpp (Maltose transport via ABC system (periplasm))
    # subsys: nothing
    # lb: 0.0, ub: 1000.0
    # (-1.0) M_atp_c + (-1.0) M_h2o_c + (-1.0) M_malt_p ==> (1.0) M_adp_c + (1.0) M_h_c + (1.0) M_malt_c + (1.0) M_pi_c
    colid = "MALTabcpp"
    @assert !hascolid(lep0, colid)
    @assert hasrowid(lep0, "malt_e")
    @assert !hasrowid(lep0, "malt_c")
    set_constraint!(lep0, colid;
        S = Dict(
            "atp_c" => -1, "h2o_c" => -1, "malt_e" => -1, 
            "adp_c" => 1, "h_c" => 1, "malt_c" => 1, "pi_c" => 1
        ),
        lb = -0.0, ub = 1000.0,
    )
    println(colid, ": ", col_str(lep0, colid))

    # NOTE
    # rxn[1725]: R_MALTptspp (Maltose transport via PEP:Pyr PTS (periplasm))
    # subsys: nothing
    # lb: 0.0, ub: 1000.0
    # (-1.0) M_malt_p + (-1.0) M_pep_c ==> (1.0) M_malt6p_c + (1.0) M_pyr_c
    # NOTE: M_malt6p_c is a dead end at iJO1366, so, using MALTabcpp

    # rxn[318]: R_AMALT1 (Amylomaltase (maltotriose))
    # subsys: nothing
    # lb: 0.0, ub: 1000.0
    # (-1.0) M_malt_c + (-1.0) M_malttr_c ==> (1.0) M_glc__D_c + (1.0) M_maltttr_c
    # I will ignore M_maltttr_c
    colid = "AMALT1"
    @assert !hascolid(lep0, colid)
    @assert hasrowid(lep0, "glc__D_e")
    @assert hasrowid(lep0, "malt_c")
    set_constraint!(lep0, colid;
        S = Dict(
            "malt_c" => -1,
            "glc__D_e" => 2,
        ),
        lb = -0.0, ub = 1000.0,
    )
    println(colid, ": ", col_str(lep0, colid))

    # -------------------------------------------
    # Glycerol
    # -------------------------------------------
    println()
    println("-"^40)
    println("Glycerol")
    println("-"^40)

    # rxn[958]: R_EX_glyc_e (Glycerol exchange)
    # subsys: nothing
    # lb: 0.0, ub: 1000.0
    # (-1.0) M_glyc_e ==> 
    colid = "EX_glyc_e"
    @assert !hascolid(lep0, colid)
    @assert !hasrowid(lep0, "glyc_e")
    set_constraint!(lep0, colid;
        S = Dict(
            "glyc_e" => -1,
        ),
        lb = -0.0, ub = 1000.0,
    )
    println(colid, ": ", col_str(lep0, colid))

    # rxn[1406]: R_GLYCtex (Glycerol transport via diffusion (extracellular to periplasm))
    # subsys: nothing
    # lb: -1000.0, ub: 1000.0
    # (-1.0) M_glyc_e <==> (1.0) M_glyc_p
    # Not apply

    # rxn[1407]: R_GLYCtpp (Glycerol transport via channel (periplasm))
    # subsys: nothing
    # lb: -1000.0, ub: 1000.0
    # (-1.0) M_glyc_c <==> (1.0) M_glyc_p
    colid = "GLYCtpp"
    @assert !hascolid(lep0, colid)
    @assert hasrowid(lep0, "glyc_e")
    @assert !hasrowid(lep0, "glyc_c")
    set_constraint!(lep0, colid;
        S = Dict(
            "glyc_e" => -1,
            "glyc_c" => 1,
        ),
        lb = -1000.0, ub = 1000.0,
    )
    println(colid, ": ", col_str(lep0, colid))

    # rxn[1408]: R_GLYK (Glycerol kinase)
    # subsys: nothing
    # lb: 0.0, ub: 1000.0
    # (-1.0) M_atp_c + (-1.0) M_glyc_c ==> (1.0) M_adp_c + (1.0) M_glyc3p_c + (1.0) M_h_c
    colid = "GLYK"
    @assert !hascolid(lep0, colid)
    @assert hasrowid(lep0, "glyc_c")
    @assert !hasrowid(lep0, "glyc3p_c")
    set_constraint!(lep0, colid;
        S = Dict(
            "atp_c" => -1, "glyc_c" => -1,
            "adp_c" => 1, "glyc3p_c" => 1, "h_c" => 1
        ),
        lb = 0.0, ub = 1000.0,
    )
    println(colid, ": ", col_str(lep0, colid))

    # rxn[1267]: R_G3PD2 (Glycerol-3-phosphate dehydrogenase (NADP))
    # subsys: nothing
    # lb: -1000.0, ub: 1000.0
    # (-1.0) M_glyc3p_c + (-1.0) M_nadp_c <==> (1.0) M_dhap_c + (1.0) M_h_c + (1.0) M_nadph_c
    colid = "G3PD2"
    @assert !hascolid(lep0, colid)
    @assert hasrowid(lep0, "glyc3p_c")
    @assert hasrowid(lep0, "dhap_c")
    @assert hasrowid(lep0, "nadp_c")
    @assert hasrowid(lep0, "nadph_c")
    set_constraint!(lep0, colid;
        S = Dict(
            "glyc3p_c" => -1, "nadp_c" => -1,
            "dhap_c" => 1, "h_c" => 1, "nadph_c" => 1
        ),
        lb = -1000.0, ub = 1000.0,
    )
    println(colid, ": ", col_str(lep0, colid))

    # -------------------------------------------
    # FBA TEST
    # -------------------------------------------
    println()
    println("="^40)
    println("FBA TEST")
    println()
    
    # EX_gal_e
    lep = emptyless_model(lep0)
    biom_id = "BIOMASS_Ecoli_core_w_GAM"
    linear_weights!(lep, biom_id, 1.0)
    lb!(lep, "EX_glc__D_e", 0.0)
    
    for nut_id in [
            "EX_glc__D_e", 
            "EX_lac__D_e", 
            "EX_malt_e", 
            "EX_gal_e", 
            "EX_glyc_e",
            "EX_ac_e"
        ]
        println("-"^40)
        println(nut_id)

        lb!(lep, nut_id, -10.0)
        opm = fba(lep, Clp.Optimizer)
        println(biom_id, ": ", solution(opm, biom_id))
        # println(nut_id, ": ", solution(opm, nut_id))
        @assert solution(opm, biom_id) > 1e-2
        lb!(lep, nut_id, 0.0)
    end    

    # -------------------------------------------
    # EXP
    # -------------------------------------------
    println()
    println("="^40)
    println("EXPERIMENTS")
    println()
    
    # -------------------------------------------
    println("-"^40)
    println("Glc regime")
    gf = 0.18
    lb!(lep, "EX_glc__D_e", -0.9 * 60 * gf)
    opm = fba(lep, Clp.Optimizer)
    println(biom_id, ": ", solution(opm, biom_id))
    lb!(lep, "EX_glc__D_e", 0.0)
    
    # -------------------------------------------
    println("-"^40)
    println("Lac-Mal-Gal regime")
    gf = 0.18
    lb!(lep, "EX_lac__D_e", -1.0 * 60 * gf)
    lb!(lep, "EX_malt_e", -0.1 * 60 * gf)
    lb!(lep, "EX_gal_e",-0.2 * 60 * gf)
    opm = fba(lep, Clp.Optimizer)
    println(biom_id, ": ", solution(opm, biom_id))
    lb!(lep, "EX_lac__D_e", 0.0)
    lb!(lep, "EX_malt_e", 0.0)
    lb!(lep, "EX_gal_e", 0.0)

    # -------------------------------------------
    println("-"^40)
    println("Glyc-Ac regime")
    gf = 0.18
    lb!(lep, "EX_glyc_e", -0.1 * 60 * gf)
    lb!(lep, "EX_ac_e", -0.6 * 60 * gf)
    opm = fba(lep, Clp.Optimizer)
    println(biom_id, ": ", solution(opm, biom_id))
    lb!(lep, "EX_glyc_e", 0.0)
    lb!(lep, "EX_ac_e", 0.0)
    
    
    # # # lb!(lep, "EX_gal_e", -0.2 * 60 * gf)
    # # # lb!(lep, "EX_glyc_e", -0.6 * 60 * gf)
    # lb!(lep, "EX_ac_e", -1.5 * 60 * gf)
    # ub!(lep, "EX_ac_e", 2.0 * 60 * gf)


end







## ------------------------------------------------------------
## ------------------------------------------------------------
## ------------------------------------------------------------
## ------------------------------------------------------------
## ------------------------------------------------------------
let
    # net
    global net0 = pull_net("iJO1366")
    
    #  open complex medium
    # Bounds from Beg et al, 2007, fig 3
    # We just need the maximum rates. 
    # Original in (mmol/ min g)
    # 1 [mmol/ min g] * 60 = 1 [mmol/ h g]
    # adjustment growth factor (I aim to adjust the maximum growth to match glucose only regime (Fig2, a))
    gf = 0.18 
    lb!(net0, "R_EX_glc__D_e", -0.9 * 60 * gf)
    # lb!(net0, "R_EX_lac__L_e", -1.0 * 60 * gf)
    # lb!(net0, "R_EX_malt_e", -0.1 * 60 * gf)
    # lb!(net0, "R_EX_gal_e", -0.2 * 60 * gf)
    # lb!(net0, "R_EX_glyc_e", -0.6 * 60 * gf)
    # lb!(net0, "R_EX_ac_e", -1.5 * 60 * gf)
    # ub!(net0, "R_EX_ac_e", 2.0 * 60 * gf)

    opm = fba(net0, Clp.Optimizer)
    rxnid = "BIOMASS_Ecoli_core_w_GAM"
    println(rxnid, ": ", solution(opm, rxnid))

end

## ------------------------------------------------------------
let
    for id in reactions(net0)
        contains(id, "EX") || continue
        summary(net0, id)
    end
end

## ------------------------------------------------------------
let
    # net
    global net0 = pull_net("ecoli_core")
    
    #  open complex medium
    # Bounds from Beg et al, 2007, fig 3
    # We just need the maximum rates. 
    # Original in (mmol/ min g)
    # 1 [mmol/ min g] * 60 = 1 [mmol/ h g]
    # adjustment growth factor (I aim to adjust the maximum growth to match glucose only regime (Fig2, a))
    gf = 0.18 
    lb!(net0, "EX_glc__D_e", -0.9 * 60 * gf)
    lb!(net0, "EX_lac__D_e", -1.0 * 60 * gf)
    # # lb!(net0, "EX_malt_e", -0.1 * 60 * gf)
    # # lb!(net0, "EX_gal_e", -0.2 * 60 * gf)
    # # lb!(net0, "EX_glyc_e", -0.6 * 60 * gf)
    # lb!(net0, "EX_ac_e", -1.5 * 60 * gf)
    # ub!(net0, "EX_ac_e", 2.0 * 60 * gf)

    opm = fba(net0, Clp.Optimizer)
    @show solution(opm, "BIOMASS_Ecoli_core_w_GAM")
    @show solution(opm, extras(net0, "EX_GLC"))

end