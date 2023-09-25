## ------------------------------------------------------------
@time begin
    using Dates
    using Random
    using MetXGEMs
    using MetXBase
    using MetXOptim
    using MetXNetHub
    using Base.Threads
    using ProgressMeter
    using NutrientLimitedGEMs
    # using CairoMakie
    using Plots
    using ProjFlows
end

## ------------------------------------------------------------
let
    N = 500
    W, H = 100.0, 100.0
    xs = ones(N) .* W/2
    ys = ones(N) .* H/2
    v = 1.0

    fns = String[]
    for ploti in 1:10

        # step!
        ploti != 1 && for ii in 1:100
            for xi in eachindex(xs)
                dx = v * rand([-1,1]) * rand()
                xs[xi] = clamp(xs[xi] + dx, 0.0, W)
            end
            for yi in eachindex(ys)
                dy = v * rand([-1,1]) * rand()
                ys[yi] = clamp(ys[yi] + dy, 0.0, W)
            end
        end

        # scatter
        ps = Plots.Plot[]
        p = plot(; xlabel = "x pos", ylabel = "y pos")
        scatter!(p, xs, ys; 
            color = :black, alpha = 0.5, 
            xlim = (0.0, W), xaxis = nothing,
            ylim = (0.0, H), yaxis = nothing,
            size = (400, 400),
            label = "", 
            # aspectratio = 1.0
        )
        push!(ps, p)

        # histogram y
        p = plot(; xlabel = "", ylabel = "")
        histogram!(p, ys; 
            label = "", color = :black, 
            bins = range(0.0, H; length = 25), 
            xaxis = nothing,
            ylim = (0.0, H), yaxis = nothing,
            orientation = :horizontal, 
            size = (400, 200),
            # aspectratio = 1.0
        )
        push!(ps, p)

        # histogram x
        p = plot(; xlabel = "", ylabel = "")
        histogram!(p, xs; 
            label = "", color = :black, 
            bins = range(0.0, W; length = 25), 
            xlim = (0.0, W), xaxis = nothing,
            yaxis = nothing,
            size = (200, 400),
            # aspectratio = 1.0
        )
        push!(ps, p)

        fn = sfig(ps, [@__DIR__, "figs"], "diff", ploti, ".png")
        push!(fns, fn)
    end

    sgif(fns, [@__DIR__], "diff.gif"; fps = 1.0) |> println

end

## ------------------------------------------------------------
## ------------------------------------------------------------
# ------------------------------------------------------------
include("1_setup_sim.jl")
include("1.1_utils.jl")
 
## ------------------------------------------------------------
# TODO: sync ECOLI_CORE and iJO1366 code
# DONE: start running the trimming script
# DONE: Finish the review
# DONE: Move to cluster SingleCell_Lactate_2023
let
    files = readdir(procdir(PROJ, [SIMVER]); join = true)
    for fn in files
        startswith(basename(fn), "obj_reg") || continue
        _, obj_reg = ldat(fn)
        do_save = false
        for obj in obj_reg
            haskey(obj, "strip.koset") || continue
            if isempty(obj["strip.koset"])
                obj["strip.koset"] = Int16[]
                do_save = true
            end
        end
        do_save || continue
        println("[", getpid(), ".", threadid(), "] ", fn)
        sdat(obj_reg, fn)
    end
end