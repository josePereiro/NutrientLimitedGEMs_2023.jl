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
    using NutrientLimitedGEMs_2023
end

# ------------------------------------------------------------
include("1_setup_sim.jl")
include("1.1_utils.jl")

# ------------------------------------------------------------
# Prepare network
@tempcontext ["KOMA" => v"0.1.0"] let

    # Hi
    _log("HELLO")

    # globals
    glob_db = query(["ROOT", "GLOBALS"])
    LP_SOLVER = glob_db["LP_SOLVER"]
    NET_ID = glob_db["NET_ID"]
    NTHREADS = glob_db["NTHREADS"]
    
    # xlep
    xlep_db = query(["ROOT", "XLEP"])
    elep0 = xlep_db["elep0"][]
    lep0 = lepmodel(elep0)
    lb0, ub0 = bounds(lep0, :)
    M, N = size(lep0)
    obj_id = extras(lep0, "BIOM")
    obj_idx = colindex(lep0, obj_id)
    obj_val_th = 0.01 # KO definition
    
    target_rxn0s = colids(lep0, elep0.idxi)
    target_rxn0is = Int16.(colindex(lep0, target_rxn0s))
    koma_hashs = UInt64[]
    # koma_idxs => status
    obj_reg = Dict{String, Any}[]
    _sync_state!(koma_hashs, obj_reg)
    
    sync_frec = 5000 # iters
    log_frec = 30.0 # seconds
    roll_count = 0
    batch_size = 5 # kos per roll
    effitiency = 1.0
    effitiency_th = 0.5 # stop if "effitiency < effitiency_th"
    
    # locks
    lkpath = procdir(PROJ, [SIMVER], "koma_sync.lk")
    thlk = ReentrantLock()
    proclk = SimpleLockFile(lkpath)
    proclk_ops = (;time_out = 10.0, recheck_time=0.1, retry_time=0.5, valid_time = 160.0, force = true)

    # koma
    @threads :static for _ in 1:NTHREADS
        
        opt_time = 0.0
        tot_time = 0.0
        init_time = time()
        last_log = 0.0
        
        th = threadid()
        th_opm = FBAOpModel(lep0, LP_SOLVER)
        set_linear_obj!(th_opm, obj_idx, MAX_SENSE)
        
        for ko in 1:Int(1e7)
            
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
                        opt_time += @elapsed optimize!(th_opm)
                        obj_val = objective_value(th_opm)
                        if obj_val > obj_val_th
                            status = :FEASIBLE
                            sort!(koset) # FEASIBLE are sorted (reduce duplication)
                        else
                            status = :UNFEASIBLE
                            _break = true;
                        end
                    end
                catch e
                    status = :ERROR
                    # _log("ERROR"; err = err_str(e; max_len = 400))
                    _break = true;
                end 
                _break && break # for toss
            end # for toss

            # up new koma
            _break = false
            lock(thlk) do
                roll_count += 1
                if status != :REVISED
                    obj = Dict{String, Any}(
                        "koset" => koset, 
                        "koma_status" => status, 
                        "koma_ver" => context("KOMA")
                    )
                    push!(obj_reg, obj)
                    # insert hash
                    i = searchsortedfirst(koma_hashs, koset_hash)
                    insert!(koma_hashs, i, koset_hash)
                end
                
                tot_time = time() - init_time
                effitiency = length(obj_reg) / roll_count

                # sync state
                if length(obj_reg) > sync_frec
                    
                    lock(proclk; proclk_ops...) do
                        _sync_state!(koma_hashs, obj_reg)
                    end
                    
                    # empty to save memmory
                    empty!(obj_reg)
                    roll_count = 0
                    GC.gc()
                end

                # LOG
                if time() - last_log > log_frec
                    lock(proclk; proclk_ops...) do
                        _log("INFO" ;
                            effitiency = effitiency,
                            opt_reltime = opt_time / tot_time,
                            roll_count = roll_count,
                            koma_hashs_len = length(koma_hashs),
                            koma_reg_len = length(obj_reg),
                            koset_len = length(koset),
                            obj_val = obj_val,
                            batch_size = batch_size,
                            koma_hashs_size = Base.summarysize(koma_hashs),
                            koma_reg_size = Base.summarysize(obj_reg),
                            status = status, 
                        )
                    end
                    last_log = time()
                end

                # if roll_count > (0.1 * sync_frec) && effitiency < effitiency_th
                #     _break = true;
                # end
            end # lock(lk) do
            _break && break # for ko
        end # for ko
    end # for _

    # save state
    lock(proclk; proclk_ops...) do
        _sync_state!(koma_hashs, obj_reg)
        _log("FINISHED")
    end

end # @tempcontext