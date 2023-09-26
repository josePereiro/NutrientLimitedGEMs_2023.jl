## ------------------------------------------------------------
# Create first a identity histogram and then a property (ej: sp) histogram.
# finally, normalize the property histogram with the weights of the identity histogram...

## ------------------------------------------------------------
let
    len = Float64[]
    col = Float64[]
    for li in keys(_obj["core_feasets"])
        val = _obj["core_feasets"][li]["core_nut_sp.EX_lac__D_e.obj_m"]
        push!(len, li)
        push!(col, val)
    end
    sperm = sortperm(len)
    lines(len[sperm], col[sperm])
end

## ------------------------------------------------------------