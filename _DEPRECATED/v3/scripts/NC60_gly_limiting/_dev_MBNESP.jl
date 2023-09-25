## ------------------------------------------------------------------
@time begin
    using MetXBase
    using MetXOptim
    using MetXNetHub
    using MetXCultureHub
    using Gurobi
    using ContextDBs
    using ExtractMacro
    using Plots
    using Serialization
end

# ------------------------------------------------------------------
include("1_setup_sim.jl")


## ------------------------------------------------------------------
## ------------------------------------------------------------------
# MAX BIOMASS NON ESSENTIALS SHADOW PRICE
function MBNESP!(opm::OpModel, biomid, essentials::Vector, targets::Vector, delta_factor)
    @assert 0.0 < delta_factor < 1.0

    # max/fix biom
    set_linear_obj!(opm, biomid, MAX_SENSE)
    optimize!(opm)
    max_biom = solution(opm, biomid)
    _biom_lb_bk, _biom_ub_bk = bounds(opm, biomid)
    bounds!(opm, biomid, max_biom, max_biom)

    # min/fix essentials intakes (assumes intakes are negative)
    set_linear_obj!(opm, essentials, MAX_SENSE) 
    optimize!(opm)
    lb_ = solution(opm, essentials)
    @assert all(lb_ .< 0)
    lb!(opm, essentials, lb_)

    # restore biom
    bounds!(opm, biomid, _biom_lb_bk, _biom_ub_bk)

    # shadow prices
    shadow_price = Float64[]
    for target in targets
        # cache target lb
        _target_lb_bk = lb(opm, target)

        # max target (assumes intakes are negative)
        set_linear_obj!(opm, target, MIN_SENSE)
        optimize!(opm)
        
        # biom0
        _target_lb0 = solution(opm, target)
        set_linear_obj!(opm, biomid, MAX_SENSE)
        optimize!(opm)
        biom0 = solution(opm, biomid)
        
        # biom1
        _target_lb1 = _target_lb0 * (1.0 - delta_factor)
        lb!(opm, target, _target_lb1)
        set_linear_obj!(opm, biomid, MAX_SENSE)
        optimize!(opm)
        biom1 = solution(opm, biomid)

        sp = (biom1 - biom0) / (_target_lb1 - _target_lb0)
        push!(shadow_price, sp)

        # restore target lb
        lb!(opm, target, _target_lb_bk)

    end

    return shadow_price
end

# ------------------------------------------------------------------
let
    # base model
    net0 = pull_net("SysBioChalmers_Human_GEM")
    
    # base model
    opm = FBAOpModel(net0, Gurobi.Optimizer)
    
    # MBNESP
    essentials = ["MAR09035", "MAR09036", "MAR09038", "MAR09039", "MAR09040", "MAR09041", "MAR09042", "MAR09043", "MAR09044", "MAR09045", "MAR09046", "MAR09048", "MAR09072", "MAR09076", "MAR09109", "MAR09143", "MAR09144", "MAR09145", "MAR09146", "MAR09151", "MAR09153", "MAR09159", "MAR09167", "MAR09269", "MAR09404"]
    non_essentials = ["MAR09034", "MAR09047", "MAR09061", "MAR09062", "MAR09063", "MAR09064", "MAR09065", "MAR09066", "MAR09067", "MAR09068", "MAR09069", "MAR09070", "MAR09071", "MAR09074", "MAR09083", "MAR09358", "MAR09361", "MAR09378", "MAR09423"]
    rx_biom = "MAR13082"
    ex_gly = "MAR09067"
    targets = essentials[1:3]
    @time prices = MBNESP!(opm, rx_biom, essentials, targets, 0.1)

    for i in eachindex(targets)
        abs(prices[i]) > 1e-5 || continue

        println("."^30)
        @show non_essentials[i]
        @show prices[i]
    end

end


## ------------------------------------------------------------------
## ------------------------------------------------------------------
## ------------------------------------------------------------------
## ------------------------------------------------------------------
## ------------------------------------------------------------------
## ------------------------------------------------------------------
let
    summary.([net0], net0.rxns[net0.subSystems .== "Isolated"]) .|> println
    return nothing
end

## ------------------------------------------------------------------
# DOING: Define ko target reactions, from subsystems
let
    target_protected_subsys = [
        "Exchange/demand reactions",
        "Artificial reactions",
    ]
end

## ------------------------------------------------------------------
let
    target_exception_subsys = [
        "Glycolysis / Gluconeogenesis",
        "Starch and sucrose metabolism",
        "Galactose metabolism",
        "Fructose and mannose metabolism",
        "Pentose and glucuronate interconversions",
        "Pyruvate metabolism",
        "Propanoate metabolism",
        "Butanoate metabolism",
        "Pentose phosphate pathway",
        "Purine metabolism",
        "Pyrimidine metabolism",
        "Nucleotide metabolism",
        "Alanine, aspartate and glutamate metabolism",
        "Arginine and proline metabolism",
        "Histidine metabolism",
        "Glycine, serine and threonine metabolism",
        "Lysine metabolism",
        "Tryptophan metabolism",
        "Tyrosine metabolism",
        "Valine, leucine, and isoleucine metabolism",
        "Phenylalanine, tyrosine and tryptophan biosynthesis",
        "Cysteine and methionine metabolism",
        "Glutathione metabolism",
        "Transport reactions",
        "Beta-alanine metabolism",
        "Metabolism of other amino acids",
        "Amino sugar and nucleotide sugar metabolism",
        "Aminoacyl-tRNA biosynthesis",
        "O-glycan metabolism",
        "N-glycan metabolism",
        "Protein assembly",
        "Protein modification",
        "Protein degradation",
        "Miscellaneous",
        "Sulfur metabolism",
        "C5-branched dibasic acid metabolism",
        "Tricarboxylic acid cycle and glyoxylate/dicarboxylate metabolism",
        "Oxidative phosphorylation",
        "ROS detoxification",
        "Fatty acid activation (cytosolic)",
        "Fatty acid activation (endoplasmic reticular)",
        "Fatty acid biosynthesis (even-chain)",
        "Fatty acid biosynthesis (odd-chain)",
        "Fatty acid biosynthesis (unsaturated)",
        "Fatty acid elongation (even-chain)",
        "Fatty acid elongation (odd-chain)",
        "Fatty acid desaturation (even-chain)",
        "Fatty acid desaturation (odd-chain)",
        "Fatty acid biosynthesis",
        "Linoleate metabolism",
        "Lipoic acid metabolism",
        "Omega-3 fatty acid metabolism",
        "Omega-6 fatty acid metabolism",
        "Arachidonic acid metabolism",
        "Leukotriene metabolism",
        "Eicosanoid metabolism",
        "Acyl-CoA hydrolysis",
        "Carnitine shuttle (cytosolic)",
        "Carnitine shuttle (mitochondrial)",
        "Carnitine shuttle (peroxisomal)",
        "Carnitine shuttle (endoplasmic reticular)",
        "Acylglycerides metabolism",
        "Beta oxidation of even-chain fatty acids (peroxisomal)",
        "Beta oxidation of odd-chain fatty acids (peroxisomal)",
        "Beta oxidation of unsaturated fatty acids (n-9) (peroxisomal)",
        "Beta oxidation of phytanic acid (peroxisomal)",
        "Beta oxidation of di-unsaturated fatty acids (n-6) (peroxisomal)",
        "Beta oxidation of even-chain fatty acids (mitochondrial)",
        "Beta oxidation of odd-chain fatty acids (mitochondrial)",
        "Beta oxidation of unsaturated fatty acids (n-7) (mitochondrial)",
        "Beta oxidation of unsaturated fatty acids (n-7) (peroxisomal)",
        "Beta oxidation of unsaturated fatty acids (n-9) (mitochondrial)",
        "Beta oxidation of di-unsaturated fatty acids (n-6) (mitochondrial)",
        "Beta oxidation of poly-unsaturated fatty acids (mitochondrial)",
        "Beta oxidation of branched-chain fatty acids (mitochondrial)",
        "Terpenoid backbone biosynthesis",
        "Steroid metabolism",
        "Androgen metabolism",
        "Cholesterol biosynthesis 1 (Bloch pathway)",
        "Cholesterol biosynthesis 2",
        "Cholesterol biosynthesis 3 (Kandustch-Russell pathway)",
        "Cholesterol metabolism",
        "Estrogen metabolism",
        "Formation and hydrolysis of cholesterol esters",
        "Sphingolipid metabolism",
        "Glycerolipid metabolism",
        "Glycerophospholipid metabolism",
        "Glycosphingolipid biosynthesis-ganglio series",
        "Glycosphingolipid biosynthesis-globo series",
        "Glycosphingolipid biosynthesis-lacto and neolacto series",
        "Glycosphingolipid metabolism",
        "Glycosylphosphatidylinositol (GPI)-anchor biosynthesis",
        "Glucocorticoid biosynthesis",
        "Prostaglandin biosynthesis",
        "Ether lipid metabolism",
        "Isolated",
        "Chondroitin / heparan sulfate biosynthesis",
        "Chondroitin sulfate degradation",
        "Heparan sulfate degradation",
        "Keratan sulfate biosynthesis",
        "Keratan sulfate degradation",
        "Bile acid biosynthesis",
        "Bile acid recycling",
        "Blood group biosynthesis",
        "Nicotinate and nicotinamide metabolism",
        "Pantothenate and CoA metabolism",
        "Inositol phosphate metabolism",
        "Folate metabolism",
        "Biotin metabolism",
        "Biopterin metabolism",
        "Ascorbate and aldarate metabolism",
        "Porphyrin metabolism",
        "Retinol metabolism",
        "Riboflavin metabolism",
        "Serotonin and melatonin biosynthesis",
        "Thiamine metabolism",
        "Vitamin B12 metabolism",
        "Vitamin B6 metabolism",
        "Vitamin D metabolism",
        "Vitamin E metabolism",
        "Xenobiotics metabolism",
        "Pool reactions",
        "Exchange/demand reactions",
        "Artificial reactions",
        "Ubiquinone synthesis",
        "Fatty acid oxidation",
        "Alkaloids biosynthesis",
        "Heme synthesis;Porphyrin metabolism",
        "Vitamin C metabolism",
        "Phosphatidylinositol phosphate metabolism",
        "Heme degradation",
        "Vitamin A metabolism",
        "Phenylalanine metabolism",
        "Heme synthesis",
        "Urea cycle",
        "Dietary fiber binding",
        "Vitamin B2 metabolism",
        "Triacylglycerol synthesis",
        "Hippurate metabolism",
        "Peptide metabolism",
        "Drug metabolism",
        "Lysine degradation",
        "Fatty acid metabolism",
        "octane oxidation",
        "	Arachidonic acid metabolism",
        "Toluene degradation",
        "Insect hormone biosynthesis",
        "Glyoxylate and dicarboxylate metabolism",
        "Fatty acid degradation",
        "Pentose Phosphate Pathway"
    ]
end



