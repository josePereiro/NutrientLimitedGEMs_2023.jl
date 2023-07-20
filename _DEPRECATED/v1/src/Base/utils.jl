## ------------------------------------------------------------------
function _deletefirst!(f::Function, col::Vector)
    idx = findfirst(f, col)
    isnothing(idx) && return idx
    deleteat!(col, idx)
    return idx
end

## ------------------------------------------------------------------
nothing