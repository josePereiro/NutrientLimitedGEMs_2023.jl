# From Chemostat_Rath2017=f511193e-9926-11ea-07ad-fb434f320181 package
# << models created by contextualizing Human1 using data from Hart et al., “High-Resolution CRISPR Screens Reveal Fitness Genes and Genotype-Specific Cancer Liabilities.” >>

const HUMAN_BIOMASS_IDER = "biomass_human"
const HUMAN_ATPM_IDER = "HMR_3964"
const PROT_POOL_EXCHANGE = "prot_pool_exchange"
const HUMAN_GLC_EX_IDER = "EX_HMR_9034"
const HUMAN_GLN_EX_IDER = "EX_HMR_9063"
const HUMAN_GAL_EX_IDER = "EX_HMR_9140"
const HUMAN_O2_EX_IDER = "EX_HMR_9048"
const HUMAN_LAC_EX_IDER = "EX_HMR_9135"
const HUMAN_AC_EX_IDER = "EX_HMR_9086"

const HART_GEM_PREFIX = "HartGEM"

function load_HartGEM_model(T=CoreModel; model_params...)

    model_file = rawdir(NutrientLimitedGEMs, HART_GEM_PREFIX, model_params, ".json")

    # load and prepare
    model = load_model(T, model_file)
    
    change_objective!(model, "biomass_human")
    # model.xl .= max.(model.xl, -500.0)
    # model.xu .= min.(model.xu, 500.0)

    
    return model
end

function check_modelfile(;params...) 
    mfile = datdir(NutrientLimitedGEMs, HART_GEM_PREFIX, params, ".json")
    exist = isfile(mfile)
    exist && @info("Model cached ", params...)
    return exist
end