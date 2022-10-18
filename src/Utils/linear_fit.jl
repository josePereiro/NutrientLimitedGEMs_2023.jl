# addapted from https://github.com/pjabardo/CurveFit.jl/blob/master/src/linfit.jl
function linear_fit(x, y)

    sx = sum(x)
    sy = sum(y)

    m = length(x)

    sx2 = zero(sx)
    sy2 = zero(sy)
    sxy = zero(sx)

    for i = 1:m
        sx2 += x[i]*x[i]
        sy2 += y[i]*y[i]
        sxy += x[i]*y[i]
    end

    a0 = (sx2*sy - sxy*sx) / ( m*sx2 - sx*sx )
    a1 = (m*sxy - sx*sy) / (m*sx2 - sx*sx)

    # serr = sum((a0 .+ a1 .* x .- y).^2) / m
    serr = sum(abs.(a0 .+ a1 .* x .- y) ./ sqrt(a1^2 + one(a1))) / m
    x0, x1 = extrema(x)
    serr = serr / ((x0 - x1)^2 + (a1 * (x0 - x1))^2)
    return (a0, a1, serr)
end

function find_lines(x, y; rtol = 1e-3, step = 2)
    lines = OrderedDict()
    len = length(x)
    atol = max(rtol, (maximum(y) - minimum(y)) * rtol)
    i0 = 1
    i = i0 + step
    while i <= len
        c, m, err = linear_fit(x[i0:i], y[i0:i])
        # @info("Doing", i0, i, err, atol, m)
        if err <= atol
            lines[i0] = (i0:i, c, m, err)
            i += step 
        else
            # new line
            i0 = i - step
            i += 1
        end
    end
    return collect(Tuple{UnitRange{Int64}, Float64, Float64, Float64}, values(lines))
end

function approx_slopes(xs, ys; n = 5)

    i0 = firstindex(xs)
    i1 = lastindex(xs)

    ms = Float64[]
    errs = Float64[]
    for i in eachindex(xs)
        r = max(i0, i - n):min(i1, i + n)
        _, m, err = linear_fit(xs[r], ys[r])
        push!(ms, m)
        push!(errs, err)
    end

    return ms, errs
end