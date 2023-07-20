@time begin
    using LinearAlgebra
    using Plots
end

# TODO: Move from here

## ------------------------------------
let
    # Generator
    g = (x) -> sin(2*π*x)
    ϵ = () -> 2e-1 * randn()
    f = (x) -> g(x) + ϵ()

    # Training Set
    N = 2500
    xs = range(0.0, 1.0; length = N)
    ts = f.(xs)

    # Polinomial approx
    M = 5 # order + 1
    λ = 0.0 # regularization

    S = zeros(M, M)
    for i in 1:M
        for j in 1:M
            for n in 1:N
                S[j, i] += xs[n]^(i - 1) * xs[n]^(j - 1)
            end
            j == 1 && continue
            S[j, i] += λ
        end
    end

    b = zeros(M)
    for j in 1:M
        for n in 1:N
            b[j] += ts[n] * xs[n]^(j-1)
        end
    end
    
    w = S\b

    f̂ = (x) -> sum(w[i]*x^(i-1) for i in 1:M)


    # Ploting
    ps = Plots.Plot[]

    p = plot(; title = "training set", xlabel = "x", ylabel = "y")
    scatter!(p, xs, ts; label = "real + noise", c = :black, alpha = 0.4, msc = :auto, m = 3)
    plot!(p, g, 0.0, 1.0; label = "real", lw = 3, c = :blue)
    plot!(p, f̂, 0.0, 1.0; label = "approx", lw = 3, c = :red)
    push!(ps, p)
    
    p = bar(w; label = "", xlabel = "j", ylabel = "wj")
    push!(ps, p)

    plot(ps...)

end