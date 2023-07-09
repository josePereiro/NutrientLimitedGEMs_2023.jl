@time begin
    using CSV
    using DataFrames
end

## ------------------------------------------------------------
let
    fn = "/Users/Pereiro/.julia/dev/NutrientLimitedGEMs/data/raw/jain_et_al_2012/CORE_data.tsv"
    df = CSV.read(fn, DataFrame)
end

## ------------------------------------------------------------
# inefficient
