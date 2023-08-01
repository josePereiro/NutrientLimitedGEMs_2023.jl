## ------------------------------------------------------------
@time begin
    using Random
    using MetXGEMs
    using MetXBase
    using MetXOptim
    using MetXNetHub
    using Base.Threads
    using ProgressMeter
    using NutrientLimitedGEMs
end

# ------------------------------------------------------------
include("1_setup_sim.jl")

## ------------------------------------------------------------
# Prepare network
@tempcontext ["KOMA" => v"0.1.0"] let

    # globals
    glob_db = query(["ROOT", "GLOBALS"])
    LP_SOLVER = glob_db["LP_SOLVER"]
    NET_ID = glob_db["NET_ID"]
    
    # xlep
    xlep_db = query(["ROOT", "XLEP"])
    global elep0 = xlep_db["elep0"][]
    global lep0 = lepmodel(elep0)
    lb0, ub0 = bounds(lep0, :)
    M, N = size(lep0)
    obj_id = extras(lep0, "BIOM")
    obj_idx = colindex(lep0, obj_id)
    obj_val_th = 0.01 # KO definition
    
    target_rxn0s = colids(lep0, elep0.idxi)
    target_rxn0is = Int16.(colindex(lep0, target_rxn0s))
    # koma_idxs => status
    _, _koma_hashs = lprocdat(PROJ, [SIMVER], "koma_hashs", ".jls") do 
        UInt64[]
    end
    global koma_hashs = _koma_hashs
    global koma_reg = Dict{Vector{Int16}, Symbol}()
    save_frec = 5000
    stop_count = 0
    batch_size = 5
    effitiency = 1.0
    effitiency_th = 0.5 # stop if "effitiency < effitiency_th"
    opt_time = 0.0
    lk = ReentrantLock()
    prog = ProgressUnknown(; dt = 1.0, desc="Progress: ", showspeed = true)


    # koma
    for _ in 1:1
        
        th = threadid()
        th_opm = FBAOpModel(lep0, LP_SOLVER)
        set_linear_obj!(th_opm, obj_idx, MAX_SENSE)
        
        for ko in 1:1500
            
            # init
            bounds!(th_opm, :, lb0, ub0) 
            th_target_rxn0is = shuffle(target_rxn0is)
            koset = Int16[]
            koset_hash = hash(0)
            
            status = :INIT
            obj_val = 0.0
            
            # coin toss
            for toss in 1:N

                # random ko
                if isempty(th_target_rxn0is) 
                    status = :EMPTY_TARGETS
                    break
                end
                for r in 1:batch_size
                    isempty(th_target_rxn0is) && break # for r
                    toko = pop!(th_target_rxn0is)
                    bounds!(th_opm, toko, 0.0, 0.0)
                    push!(koset, toko)
                end
                koset_hash = hash(koset)

                _break = false
                try
                    # Test koma
                    if insorted(koset_hash, koma_hashs)
                        status = :REVISED
                        _break = true;
                    else
                        # Test biomass
                        opt_time = @elapsed optimize!(th_opm)
                        obj_val = objective_value(th_opm)
                        if obj_val > obj_val_th
                            status = :FEASIBLE
                        else
                            status = :UNFEASIBLE
                            _break = true;
                        end
                    end
                catch e
                    status = :ERROR
                    showerror(stdout, e); println(); flush(stdout);
                    _break = true;
                end 
                _break && break # for toss
            end # for toss

            # up new koma
            _break = false
            lock(lk) do
                stop_count += 1
                if status != :REVISED
                    koma_reg[koset] = status
                    # insert hash
                    i = searchsortedfirst(koma_hashs, koset_hash)
                    insert!(koma_hashs, i, koset_hash)
                end
                effitiency = length(koma_reg) / stop_count
                
                next!(prog; showvalues = () -> [
                    (:th, th), 
                    (:opt_time, opt_time),
                    (:effitiency, effitiency),
                    (:koma_hashs_len, length(koma_hashs)),
                    (:koma_reg_len, length(koma_reg)),
                    (:batch_size, batch_size),
                    (:koma_hashs_size, Base.summarysize(koma_hashs)),
                    (:koma_reg_size, Base.summarysize(koma_reg)),
                    (:status, status), 
                    (:koset, join(koset, ", ")),
                ])

                # save state
                if length(koma_reg) > save_frec
                    sprocdat(PROJ, koma_hashs, 
                        [SIMVER], "koma_hashs", ".jls"
                    )
                    sprocdat(PROJ, koma_reg, 
                        [SIMVER], "koma_reg", (;h = hash(koma_hashs)), ".jls"
                    )
                    # empty to save memmory
                    empty!(koma_reg)
                    stop_count = 0
                    GC.gc()
                end

                if length(koma_reg) > 0.1 * save_frec && effitiency < effitiency_th 
                    _break = true;
                end
            end # lock(lk) do
            _break && break # for ko
        end # for ko
    end # for _


    # save state
    sprocdat(PROJ, koma_hashs, 
        [SIMVER], "koma_hashs", ".jls"
    )
    sprocdat(PROJ, koma_reg, 
        [SIMVER], "koma_reg", (;h = hash(koma_hashs)), ".jls"
    )

    
    finish!(prog)

end # @tempcontext