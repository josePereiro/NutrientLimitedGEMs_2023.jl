## ------------------------------------------------------------------
function _collect_ill_trajs(traj_dir, alg, ΔSth)
    illtraj_fns = String[]
    @time for fn in readdir(traj_dir; join = true)
        
        endswith(fn, ".jls") || continue
        traj = ldat(fn)
        traj["status"] == :success || continue
        epdat = get(traj, alg, nothing)
        isnothing(epdat) && continue

        # Ss & V
        Ss = epdat["Ss"]
        
        # ep_status1
        ep_statuses = epdat["ep_statuses"]
        ep_status1 = last(ep_statuses)
        all(ep_statuses .== :converged) || continue

        # traj_idxs
        traj_idxs = traj["traj_idxs"]
        isempty(traj_idxs) && continue

        Ssis = findall(Ss) do S
            iszero(S) && return false
            isnan(S) && return false
            return true
        end
        isempty(Ssis) && continue

        ΔS = last(Ss[Ssis]) - first(Ss[Ssis])
        
        if ΔS > ΔSth
            push!(illtraj_fns, fn)
            @show ΔS
        end

    end
    return illtraj_fns
end

## ------------------------------------------------------------------
function _merge_hist_pool!(hist_pool::Dict, net, rxn, ws, samples; 
        nbins = 100, margin = 0.3
    )
    rxni = rxnindex(net, rxn)
    rxn = net.rxns[rxni]

    # up histogram
    bins, hist = get!(hist_pool, rxn) do
        Δv = rxnrange(net, rxni)
        l, u = bounds(net, rxni)
        _bins = range(l - margin * Δv, u + margin * Δv; length = nbins)
        _hist = zeros(nbins)
        return (_bins, _hist)
    end
    _histogram!(bins, hist, view(samples, :, rxni), ws)
end

## ------------------------------------------------------------------
function _reduce(alg)
    
    samples_dir = joinpath(@__DIR__, "hists")

    hr_hist_pool0 = Dict()
    epmc_hist_pool0 = Dict()
    ep_join_hist_pool0 = Dict()
    ep_boxed_hist_pool0 = Dict()
    ep_tilted_hist_pool0 = Dict()

    net = nothing
    for file in readdir(samples_dir; join = true)
        endswith(file, ".jls") || continue
        basename(file) == "histograms.jls" && continue
        
        dat = ldat(file)
        dat.alg == alg || continue
        
        @show file
        
        net = dat.net
        
        for (hist_pool, hist_pool0) in [
                (dat.hr_hist_pool, hr_hist_pool0),
                (dat.epmc_hist_pool, epmc_hist_pool0),
                (dat.ep_join_hist_pool, ep_join_hist_pool0),
                (dat.ep_boxed_hist_pool, ep_boxed_hist_pool0),
                (dat.ep_tilted_hist_pool, ep_tilted_hist_pool0),
            ]

            if isempty(hist_pool0)
                merge!(hist_pool0, hist_pool)
                continue
            end

            for rxn in net.rxns
                bins0, hist0 = hist_pool0[rxn]
                bins, hist = hist_pool[rxn]
                hist0 .+= hist
                @assert all(bins .== bins0)
            end
        end
        
    end

    isnothing(net) && return

    dat = (;
        alg, net, 
        hr_hist_pool = hr_hist_pool0,
        epmc_hist_pool = epmc_hist_pool0,
        ep_join_hist_pool = ep_join_hist_pool0,
        ep_boxed_hist_pool = ep_boxed_hist_pool0,
        ep_tilted_hist_pool = ep_tilted_hist_pool0,
    )
    
    return dat

end

function _tryreduce(alg)
    try; _reduce(alg)        
        catch err; @show err
    end
end

## ------------------------------------------------------------------
# Collect Marginals
function _accumulate_marginals(net, alg::String; 
        nbins = 100,
        range_margin = 0.3,
        nsamples = Int(1e5),
        runid = hash(time()),
        save_frec = 60.0, # secs
        tout = Inf,  # secs
    )

    hr_hist_pool = Dict()
    epmc_hist_pool = Dict()
    ep_join_hist_pool = Dict()
    ep_boxed_hist_pool = Dict()
    ep_tilted_hist_pool = Dict()

    # -------------------------------
    # EP
    epm = FluxEPModelT0(net)
    config!(epm; 
        damp = 0.9,
        verbose = true,
        epsconv = 1e-6, 
        maxiter = Int(1e6)
    )
    converge!(epm)
    @show convergence_status(epm)
    
    # -------------------------------
    # SAMPLERS
    hrm = HRModel(net, Gurobi.Optimizer) # HitAndRun
    epmc = EPBoxedMvNormal(epm) # Uniform (EPBoxedMvNormal)
    epjd = EPJoinMvNormal(epm) # EPJoinMvNormal
    epbd = EPBoxedMvNormal(epm) # EPBoxedMvNormal
    etbds = [EPTiltedMvNormal(epm, rxn) for rxn in net.rxns] # EPTiltedMvNormal
    
    t0, tsave = time(), time()
    while true
        
        # -------------------------------
        # HR SAMPLES
        println("."^40); @info("HRMC")
        @time hr_ws, hr_samples = sample!(hrm, nsamples; rw = one)
        for rxn in net.rxns
            _merge_hist_pool!(hr_hist_pool, net, rxn, hr_ws, hr_samples; 
                nbins, margin = range_margin
            )
        end
        hr_samples = nothing

        # # EPMC SAMPLES
        # println("."^40); @info("EPMC")
        # @time epmc_ws, epmc_samples = sample!(epmc, 5 * nsamples; rw = inv)
        # for rxn in net.rxns
        #     _merge_hist_pool!(epmc_hist_pool, net, rxn, epmc_ws, epmc_samples; 
        #         nbins, margin = range_margin
        #     )
        # end
        # epmc_samples = nothing

        # # EP BOXED JOIN SAMPLES
        # println("."^40); @info("EP BOXED")
        # @time ep_boxed_ws, ep_boxed_samples = sample!(epbd, nsamples; rw = inv)
        # for rxn in net.rxns
        #     _merge_hist_pool!(ep_boxed_hist_pool, net, rxn, ep_boxed_ws, ep_boxed_samples; 
        #         nbins, margin = range_margin
        #     )
        # end
        # ep_boxed_samples = nothing
        
        # EP JOIN SAMPLES
        println("."^40); @info("EP JOIN")
        @time ep_join_ws, ep_join_samples = sample!(epjd, nsamples; rw = one)
        for rxn in net.rxns
            _merge_hist_pool!(ep_join_hist_pool, net, rxn, ep_join_ws, ep_join_samples; 
                nbins, margin = range_margin
            )
        end
        ep_join_samples = nothing

        # EP TILTED SAMPLES
        for (rxn, etbd) in zip(net.rxns, etbds)
            println("."^40); @info("EP TILTED", rxn)
            @time ep_tilted_ws, ep_tilted_samples = sample!(etbd, nsamples; rw = one)
            _merge_hist_pool!(ep_tilted_hist_pool, net, rxn, ep_tilted_ws, ep_tilted_samples; 
                nbins, margin = range_margin
            )
            ep_tilted_samples = nothing
        end

        # -------------------------------
        # SAVE
        abs(tsave - time()) > save_frec || continue
        tsave = time()
        
        fn = dfname([@__DIR__, "hists"], runid, ".jls")
        
        println("."^40)
        @info("Saved!!!", save_frec, fn = basename(fn))
        dat = (;
            alg, net, 
            hr_hist_pool,
            epmc_hist_pool,
            ep_join_hist_pool,
            ep_boxed_hist_pool,
            ep_tilted_hist_pool,
        )
        sdat(dat, fn)

        # -------------------------------
        # tout
        abs(t0 - time()) > tout && break

        # -------------------------------
        # REDO MODELS (Just in case)
        hrm = HRModel(net, Gurobi.Optimizer) # HitAndRun
        epmc = EPBoxedMvNormal(epm) # Uniform (EPBoxedMvNormal)
        epjd = EPJoinMvNormal(epm) # EPJoinMvNormal
        epbd = EPBoxedMvNormal(epm) # EPBoxedMvNormal
        etbds = [EPTiltedMvNormal(epm, rxn) for rxn in net.rxns] # EPTiltedMvNormal

    end
end

## ------------------------------------------------------------------
function _plot_hists(alg)
    
    println("."^40)
    @info("Reducing")
    dat = _reduce(alg)
    isnothing(dat) && return
    net = dat.net

    # EPMC
    println("."^40)
    @info("Ploting")

    ps = Plots.Plot[]
    for rxn in sort(net.rxns)
        
        @show rxn
        p = plot(;xlabel = rxn, ylabel = "~ pdf")

        x0, x1 = Inf, -Inf
        for (hist_pool, label) in [
            (dat.hr_hist_pool, "HRMC"),
            (dat.epmc_hist_pool, "EPMC"),
            (dat.ep_join_hist_pool, "EP Join"),
            # (dat.ep_boxed_hist_pool, "EP Boxed"),
            (dat.ep_tilted_hist_pool, "EP Tilted"),
        ]
            isempty(hist_pool) && continue

            bins, hist0 = hist_pool[rxn]
            hist = hist0 ./ sum(hist0)
            plot!(p, bins, hist; label, lw = 3)
            # scatter!(p, bins, hist; label, m = 3, msc = :auto)
            
            # xlim detection
            hist = (hist0 .- minimum(hist0)) 
            hist = hist ./ maximum(hist)

            highbins = bins[hist .> 1e-3]
            x0 = min(first(highbins), x0)
            x1 = max(last(highbins), x1)
        end

        # plot!(; xlim = [x0, x1])
        
        push!(ps, p)

    end

    # Save
    println("."^40)
    @info("Saving")

    fn = dfname(
        [@__DIR__, "plots"], 
        alg, "marginals", ".png"
    )
    sfig(ps, fn)
    @show fn

end