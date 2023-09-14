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