## ------------------------------------------------------------------
# Experimental data
_, (CalDat, _) = lprocdat(PROJ,
    ["Calzadilla_et_al"], "Calzadilla.bundle", ".jls"; 
    verbose = false
)
nothing

## ------------------------------------------------------------------
_Calzadilla_hek_idmap = Dict(
    "ALA" => "R_EX_ala_L_LPAREN_e_RPAREN_",
    "ARG" => "R_EX_arg_L_LPAREN_e_RPAREN_",
    "ASN" => "R_EX_asn_L_LPAREN_e_RPAREN_",
    "ASP" => "R_EX_asp_L_LPAREN_e_RPAREN_",
    "GLN" => "R_EX_gln_L_LPAREN_e_RPAREN_",
    "GLU" => "R_EX_glu_L_LPAREN_e_RPAREN_",
    "GLC" => "R_EX_glc_LPAREN_e_RPAREN_",
    "GLY" => "R_EX_gly_LPAREN_e_RPAREN_",
    "ILE" => "R_EX_ile_L_LPAREN_e_RPAREN_",
    "LAC" => "R_EX_lac_L_LPAREN_e_RPAREN_",
    "LEU" => "R_EX_leu_L_LPAREN_e_RPAREN_",
    "LYS" => "R_EX_lys_L_LPAREN_e_RPAREN_",
    "MET" => "R_EX_met_L_LPAREN_e_RPAREN_",
    "PHE" => "R_EX_phe_L_LPAREN_e_RPAREN_",
    "PRO" => "R_EX_pro_L_LPAREN_e_RPAREN_",
    "SER" => "R_EX_ser_L_LPAREN_e_RPAREN_",
    "THR" => "R_EX_thr_L_LPAREN_e_RPAREN_",
    "TRP" => "R_EX_trp_L_LPAREN_e_RPAREN_",
    "TYR" => "R_EX_tyr_L_LPAREN_e_RPAREN_",
    "VAL" => "R_EX_val_L_LPAREN_e_RPAREN_",
)
for (k,v) in _Calzadilla_hek_idmap
    _Calzadilla_hek_idmap[v] = k
end

## ------------------------------------------------------------------
function _Calzadilla_heknet!(lep0, st;
        solver = LP_SOLVER,
    )

    # contextualization
    # lb = - c D / Xv

    stdat = CalDat[st]
    D = stdat["D"]["val"]
    Xv = stdat["Xv"]["val"]

    # ["ALA", "ARG", "ASN", "ASP", "GLN", "GLU", "GLC", 
    # "GLY", "ILE", "LAC", "LEU", "LYS", "MET", "PHE", 
    # "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
    for CalId in CalDat["metids"]

        c = stdat["c$(CalId)"]["val"]
        lb = -(c * D / Xv)
        
        HekId = _Calzadilla_hek_idmap[CalId]
        # println(CalId, " lb: ", lb)
        lb!(lep0, HekId, lb)

    end
    
    # TODO: ask Calzadilla
    lb!(lep0, "R_EX_tyr_L_LPAREN_e_RPAREN_", -1.0) # To restrictive
    lb!(lep0, "R_EX_arg_L_LPAREN_e_RPAREN_", -1.0) # To restrictive
    lb!(lep0, "R_EX_his_L_LPAREN_e_RPAREN_", -1.0)
    
    lb!(lep0, "R_EX_o2_LPAREN_e_RPAREN_", -1000.0)
    lb!(lep0, "R_EX_h_LPAREN_e_RPAREN_", -1000.0)
    lb!(lep0, "R_EX_pi_LPAREN_e_RPAREN_", -1000.0)
    lb!(lep0, "R_EX_h2o_LPAREN_e_RPAREN_", -1000.0)
    lb!(lep0, "R_EX_co2_LPAREN_e_RPAREN_", -1000.0)
    lb!(lep0, "R_EX_hco3_LPAREN_e_RPAREN_", -1000.0)

    # max biom
    # opm = FBAOpModel(lep0, solver)
    # optimize!(opm)
    # biom_id = extras(lep0, "BIOM")
    # biom0 = solution(opm, biom_id)
    # @show biom0
    # @show D

    # _, lep = withcachedat(PROJ, :get!, (lep0,);
    #     verbose = true
    # ) do 
    #     box(lep0, solver; verbose = true)
    # end
    return lep0
end

function _Calzadilla_heknet(arg...; kwargs...)
    model_id = "Martinez_Monge_HEK293"
    net = pull_net(model_id)
    clampbounds!(net, -1000.0, 1000.0)
    _Calzadilla_heknet!(net, arg...; kwargs...)
end

## ------------------------------------------------------------------
# Martinez_Monge_HEK293
function _base_heknet(PROJ::Project0; 
        solver = LP_SOLVER,
    )

    model_id = "Martinez_Monge_HEK293"
    net = pull_net(model_id)
    clampbounds!(net, -1000.0, 1000.0)

    # old contextualization
    lb!(net, "R_EX_tyr_L_LPAREN_e_RPAREN_", -1.0)
    lb!(net, "R_EX_his_L_LPAREN_e_RPAREN_", -1.0)
    lb!(net, "R_EX_asn_L_LPAREN_e_RPAREN_", -1.0)
    lb!(net, "R_EX_arg_L_LPAREN_e_RPAREN_", -1.0)
    lb!(net, "R_EX_ile_L_LPAREN_e_RPAREN_", -1.0)
    lb!(net, "R_EX_phe_L_LPAREN_e_RPAREN_", -1.0)
    lb!(net, "R_EX_met_L_LPAREN_e_RPAREN_", -1.0)
    lb!(net, "R_EX_thr_L_LPAREN_e_RPAREN_", -1.0)
    lb!(net, "R_EX_lys_L_LPAREN_e_RPAREN_", -1.0)
    lb!(net, "R_EX_leu_L_LPAREN_e_RPAREN_", -1.0)
    lb!(net, "R_EX_val_L_LPAREN_e_RPAREN_", -1.0)
    lb!(net, "R_EX_trp_L_LPAREN_e_RPAREN_", -1.0)

    lb!(net, "R_EX_glc_LPAREN_e_RPAREN_", -1.0)

    lb!(net, "R_EX_glu_L_LPAREN_e_RPAREN_", -1.0)
    lb!(net, "R_EX_gln_L_LPAREN_e_RPAREN_", -1.0)
    # lb!(net, "R_EX_lac_L_LPAREN_e_RPAREN_", -1000.0)
    # lb!(net, "R_EX_nh4_LPAREN_e_RPAREN_", -1000.0)

    lb!(net, "R_EX_o2_LPAREN_e_RPAREN_", -1000.0)
    lb!(net, "R_EX_h_LPAREN_e_RPAREN_", -1000.0)
    lb!(net, "R_EX_pi_LPAREN_e_RPAREN_", -1000.0)
    lb!(net, "R_EX_h2o_LPAREN_e_RPAREN_", -1000.0)
    lb!(net, "R_EX_co2_LPAREN_e_RPAREN_", -1000.0)
    lb!(net, "R_EX_hco3_LPAREN_e_RPAREN_", -1000.0)

    # max biom
    # opm = FBAOpModel(net, solver)
    # optimize!(opm)
    # biom_id = extras(net, "BIOM")
    # biom0 = solution(opm, biom_id)
    # @show biom0
    # # @show D

    _, lep = withcachedat(PROJ, :get!, (net,);
        verbose = true
    ) do 
        box(net, solver; verbose = true)
    end
    return lep
end

## ------------------------------------------------------------------
function _foreach_sim(f::Function, PROJ, netid, simid)
    
    traj_dir = procdir(PROJ, [netid, "sims"])

    @time for fn in readdir(traj_dir; join = true)
        endswith(fn, ".jls") || continue
        _, sim = ldat(fn)
        sim["status"] == :success || continue
        haskey(sim, simid) || continue
        
        epdat = sim[simid]

        # stuff
        Fs = epdat["Fs"]
        Ss = epdat["Ss"]
        ep_statuses = epdat["ep_statuses"]

        val_idxs1 = findall(ep_statuses .== :converged)
        val_idxs2 = NutrientLimitedGEMs._find_val_idxs(Ss, Fs)
        val_idxs = intersect(val_idxs1, val_idxs2)
        isempty(val_idxs) && continue

        f(sim, epdat, val_idxs) === true && break
    end
end

## ------------------------------------------------------------------
subSyst_map = Dict(
    "Transport, mitochondrial" => "Transport",
    "Glycine, serine, alanine and threonine metabolism" => "AA metabolism",
    "Lysine metabolism" => "AA metabolism",
    "Tryptophan metabolism" => "AA metabolism",
    "Sphingolipid metabolism" => "Lipid metabolism",
    "Valine, leucine, and isoleucine metabolism" => "AA metabolism",
    "Transport, extracellular" => "Transport",
    "Fatty acid oxidation" => "Lipid metabolism",
    "Citric acid cycle" => "Citric acid cycle",
    "Nucleotide interconversion" => "Nucleotide metabolism",
    "Purine synthesis" => "AA metabolism",
    "Triacylglycerol synthesis" => "Lipid metabolism",
    "Methionine and cysteine metabolism" => "AA metabolism",
    "null" => "",
    "Glutamate metabolism" => "AA metabolism",
    "Urea cycle" => "Urea cycle",
    "R group synthesis" => "AA metabolism",
    "Pyrimidine synthesis" => "Pyrimidine synthesis",
    "Alanine and aspartate metabolism" => "AA metabolism",
    "Oxidative phosphorylation" => "Oxidative phosphorylation",
    "Cholesterol metabolism" => "Lipid metabolism",
    "Glycerophospholipid metabolism" => "Lipid metabolism",
    "Transport, endoplasmic reticular" => "Transport",
    "Transport, nuclear" => "Transport",
    "Folate metabolism" => "Folate metabolism",
    "Exchange/demand reaction" => "Transport",
    "Squalene and cholesterol synthesis" => "Lipid metabolism",
    "Glycolysis/gluconeogenesis" => "Sugar metabolism",
    "" => "",
    "Fatty acid synthesis" => "Lipid metabolism",
    "Arginine and Proline Metabolism" => "AA metabolism",
    "Pentose phosphate pathway" => "Sugar metabolism",
    "Nucleotide salvage pathway" => "Nucleotide metabolism",
    "Histidine metabolism" => "AA metabolism",
    "Miscellaneous" => "Miscellaneous",
    "Pyruvate metabolism" => "Pyruvate metabolism",
    "Inositol phosphate metabolism" => "Inositol phosphate metabolism",
    "Purine catabolism" => "Nucleotide metabolism",
    "Cysteine Metabolism" => "AA metabolism",
)

subSyst_colors = begin 
    subs = unique(collect(values(subSyst_map)))
    colors = Plots.distinguishable_colors(length(subs))
    Dict(subs .=> colors)
end


nothing