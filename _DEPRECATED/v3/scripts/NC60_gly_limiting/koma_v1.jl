## ------------------------------------------------------------
@time begin
    using CSV
    using DataFrames
    using MetXBase
    using MetXOptim
    using MetXNetHub
    using Base.Threads
    using Clp
    using Random
    using ProgressMeter
    using UnicodePlots
end

## ------------------------------------------------------------
# Generate sequences
let
    net0 = pull_net("ecoli_core")
    # net0 = pull_net("ECC2")
    # net0 = pull_net("iCHO2291")
    # net0 = pull_net("SysBioChalmers_Human_GEM")
    lep0 = lepmodel(net0)
    lb0, ub0 = bounds(lep0, :)
    M, N = size(lep0)
    obj_id = extras(lep0, "BIOM")
    obj_idx = colindex(lep0, obj_id)
    obj_val_th = 0.01
    
    target_rxn0s = colids(lep0)
    target_rxn0is = Int16.(colindex(lep0, target_rxn0s))
    global koma = Dict{Set{Int16}, Symbol}()
    stop_count = 0
    batch_size = 5
    effitiency = 1.0
    effitiency_th = 0.5 # stop if "effitiency < effitiency_th"
    lk = ReentrantLock()
    prog = ProgressUnknown(; dt = 1.0, desc="Progress: ")

    # generator
    @threads :static for _ in 1:nthreads()
    # for _ in 1:1
        
        th = threadid()
        th_opm = FBAOpModel(lep0, Clp.Optimizer)
        set_linear_obj!(th_opm, obj_idx, MAX_SENSE)
        
        for ko in 1:500
            
            # init
            bounds!(th_opm, :, lb0, ub0) 
            th_target_rxn0is = shuffle(target_rxn0is)
            koset = Set{Int16}()
            
            status = :INIT
            obj_val = 0.0
            
            # coin toss
            for toss in 1:N

                # random ko
                isempty(th_target_rxn0is) && (status = :EMPTY_TARGETS; break)
                for r in 1:rand(1:batch_size)
                    isempty(th_target_rxn0is) && break
                    toko = pop!(th_target_rxn0is)
                    bounds!(th_opm, toko, 0.0, 0.0)
                    push!(koset, toko)
                end

                try
                    # Test koma
                    if haskey(koma, koset)
                        status = koma[koset]
                        break
                    end

                    # Test biomass
                    optimize!(th_opm)
                    obj_val = objective_value(th_opm)
                    if obj_val > obj_val_th
                        status = :FEASIBLE
                    else
                        status = :UNFEASIBLE
                        break
                    end
                catch e
                    status = :ERROR
                    # showerror(stdout, e); println(); flush(stdout);
                    break
                end
            end

            lock(lk) do
                stop_count += 1
                koma[koset] = status
                effitiency = length(koma) / stop_count
                
                next!(prog; showvalues = () -> [
                    (:koma_length, length(koma)),
                    (:effitiency, effitiency),
                    (:batch_size, batch_size),
                    (:koma_size, Base.summarysize(koma)),
                    (:th, th), 
                    (:status, status), 
                    (:koset, join(koset, ", ")),
                    # (:hist, strip(string(UnicodePlots.histogram(length.(keys(koma)); ylabel = "ko count", xlabel = "frec")))), 
                ])

                if effitiency < effitiency_th 
                    break;
                end
            end

        end # for ko in 1:3

    end
    
    finish!(prog)

    
end

## ------------------------------------------------------------
let
    lens = length.(keys(koma))
    str = string(UnicodePlots.histogram(lens; xlabel = "ko count", ylabel = "frec"))
    println(str)
end

## ------------------------------------------------------------
function haspattern(set::Vector{Bool}, pattern::Vector{Bool})
    @inbounds for (ti, si) in zip(set, pattern)
        si && !ti && return false
    end
    return true
end

# ## ------------------------------------------------------------
# let
#     fn = "/Users/Pereiro/.julia/dev/NutrientLimitedGEMs_2023/data/raw/jain_et_al_2012/CORE_data.tsv"
#     df = CSV.read(fn, DataFrame)
# end
