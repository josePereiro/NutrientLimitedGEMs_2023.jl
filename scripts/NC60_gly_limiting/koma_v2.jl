## ------------------------------------------------------------
@time begin
    using CSV
    using DataFrames
    using MetXBase
    using MetXOptim
    using MetXNetHub
    using MetXGEMs
    using Base.Threads
    using Clp
    using Random
    using ProgressMeter
    using Plots
    using UnicodePlots
end

## ------------------------------------------------------------
let
    println()
    println("="^60)
    println("IO TEST")
    println("."^60)
    println()

    # setup
    test_dir = tempdir()
    mkpath(test_dir)
    
    # net0
    lep0 = MetXGEMs.toy_model() |> lepmodel
    
    for ext in [".mat", 
        # ".xml", # TODO: Check why this error (COBREXA related)
        ".json"]
        
        fn = joinpath(test_dir, string("test_net", ext))
        @show basename(fn)
        try
            rm(fn; force = true)
            # save
            save_net(lep0, fn)
            @assert isfile(fn)
            
            # load
            lep1 = load_net(fn)

            # Test LEP fields
            @assert all(colids(lep0) .== colids(lep1))
            @assert all(rowids(lep0) .== rowids(lep1))
            @assert all(cost_matrix(lep0) .== cost_matrix(lep1))
            @assert all(balance(lep0) .== balance(lep1))
            @assert all(lb(lep0) .== lb(lep1))
            @assert all(ub(lep0) .== ub(lep1))
        finally
            rm(fn; force = true)
        end

    end


    println()
end

## ------------------------------------------------------------
# Prepare network
# TODO: Do this and store it at MetXNetHub (document the source version)
let
    # netid = "ecoli_core"
    # netid = "iJR904"
    netid = "SysBioChalmers_Human_GEM"
    net0 = pull_net(netid)
    lep0 = lepmodel(net0)
    # lep1 = box(lep0, Clp.Optimizer)
    lep1 = lep0
    # @time global elep1 = EchelonLEPModel(lep1; verbose = true)
    @time global elep1 = lep1
    # @show length(elep1.idxi)
    @show size(lep1, 2)

    for ext in [".xml", ".mat"]
        fn = joinpath(@__DIR__, string("elep1-",netid, ext))
        @show fn
        save_net(elep1, fn)
    end
    nothing
end

## ------------------------------------------------------------
# Generate sequences
let
    # net0 = pull_net("toy_net")
    net0 = pull_net("ecoli_core")
    lep0 = lepmodel(net0)
    enet0 = EchelonLEPModel(lep0)
    # TODO use EchelonForm for work with independent only
    # net0 = pull_net("ECC2")
    # net0 = pull_net("iCHO2291")
    # net0 = pull_net("SysBioChalmers_Human_GEM")
    lb0, ub0 = bounds(lep0, :)
    M, N = size(lep0)
    obj_id = extras(lep0, "BIOM")
    obj_idx = colindex(lep0, obj_id)
    obj_val_th = 0.01
    

    target_rxn0s = colids(lep0)
    target_rxn0is = Set{Int16}(colindex(lep0, target_rxn0s))
    lk = ReentrantLock()
    # prog = ProgressUnknown(; dt = 1.0, desc="Progress: ")

    # generator state
    idle_counters = Dict{Int, Int}()
    global seeds = Set([Set{Int16}([ri]) for ri in target_rxn0is]) # set initial seeds
    global unfeasibles = Dict{Set{Int16}, Symbol}()

    # params
    koset_max = 3000

    # generator
    @threads for _ in 1:3*nthreads()
    # for _ in 1:1
        th = threadid()

        th_opm = FBAOpModel(lep0, Clp.Optimizer)
        set_linear_obj!(th_opm, obj_idx, MAX_SENSE)

        for s in 1:100000000
        
            tocheck = Vector{Set{Int16}}()

            lock(lk) do
                println("."^60)
                @show th

                # idle handling
                if isempty(seeds)
                    counter = get!(idle_counters, th, 0)
                    idle_counters[th] += 1
                    @info("Idle"); 
                    counter = idle_counters[th]
                    @show counter
                    sleep(0.5)
                    return
                end
                
                # compute to check
                @show length(seeds)
                seed = pop!(seeds)
                remainders = setdiff(target_rxn0is, seed)
                @show length(remainders)
                if isempty(remainders) # nothing left to grow
                    seed = nothing
                end
                tocheck = [push!(copy(seed), remainder) for remainder in remainders]
                @show length(tocheck)
                @show length(unfeasibles)
                
            end
            get(idle_counters, th, 0) > 3 && break # loop

            # ko
            for koset in tocheck
                
                # init
                bounds!(th_opm, :, lb0, ub0)
                feasible = rand() > 0.2

                # ko
                for ko in koset
                    bounds!(th_opm, ko, 0.0, 0.0)
                end
                
                try
                    # Test biomass
                    optimize!(th_opm)
                    obj_val = objective_value(th_opm)
                    if obj_val > obj_val_th
                        status = :FEASIBLE
                        if length(koset) <= koset_max
                            push!(seeds, koset) # keep going
                        end # ignore large kosets
                    else
                        unfeasibles[koset] = :UNFEASIBLE
                    end
                catch e
                    unfeasibles[koset] = :ERROR
                    # showerror(stdout, e); println(); flush(stdout);
                end

            end # for koset

        end # for it
    
    
    end # @threads
    
    # finish!(prog)

    
end

## ------------------------------------------------------------
let
    lens = length.(keys(koma))
    str = string(UnicodePlots.histogram(lens; xlabel = "ko count", ylabel = "frec"))
    println(str)
end

## ------------------------------------------------------------

## ------------------------------------------------------------
function haspattern(set::Vector{Bool}, pattern::Vector{Bool})
    @inbounds for (ti, si) in zip(set, pattern)
        si && !ti && return false
    end
    return true
end

# ## ------------------------------------------------------------
# let
#     fn = "/Users/Pereiro/.julia/dev/NutrientLimitedGEMs/data/raw/jain_et_al_2012/CORE_data.tsv"
#     df = CSV.read(fn, DataFrame)
# end
