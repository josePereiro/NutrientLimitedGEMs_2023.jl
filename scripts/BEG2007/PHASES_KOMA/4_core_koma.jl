## ------------------------------------------------------------
@time begin
    using MetXGEMs
    using MetXBase
    using MetXOptim
    using MetXNetHub
    using BlobBatches
    using Base.Threads
    using NutrientLimitedGEMs
end

# ------------------------------------------------------------
include("1_setup.jl")
include("2_utils.jl")

## ------------------------------------------------------------
# Prepare network
@tempcontext ["CORE_KOMA" => v"0.1.0"] let

    # context
    ALG_VER = context("CORE_KOMA")

    # Hi
    println("[", getpid(), ".", threadid(), "] ", "HELLO")

    # globals
    GLOB_DB = query(["ROOT", "GLOBALS"])
    LP_SOLVER = GLOB_DB["LP_SOLVER"]
    NTHREADS = GLOB_DB["NTHREADS"]
    
    # xlep
    core_xlep_db = query(["ROOT", "CORE_XLEP"])
    core_elep0 = core_xlep_db["core_elep0"][]
    core_lep0 = lepmodel(core_elep0)
    lb0, ub0 = bounds(core_lep0, :)
    M, N = size(core_lep0)
    obj_id = extras(core_lep0, "BIOM")
    obj_idx = colindex(core_lep0, obj_id)
    obj_val_th = 0.01 # KO definition
    
    target_rxn0s = colids(core_lep0, core_elep0.idxi)
    target_rxn0is = Int16.(colindex(core_lep0, target_rxn0s))
    koma_hashs = UInt64[]
    _sync_koma_hashs!(koma_hashs)
    
    max_nblobs = 500 # bobs per batch
    log_frec = 10.0 # seconds
    roll_count = 0
    batch_size = 3 # kos per roll
    effitiency = 1.0
    effitiency_th = -1 # stop if "effitiency < effitiency_th"
    DOWNREG_FACTOR = GLOB_DB["DOWNREG_FACTOR"] # stop if "effitiency < effitiency_th"

    # koma
    opt_time = 0.0
    tot_time = 0.0
    init_time = time()
    last_log = 0.0
    
    th = threadid()
    th_opm = FBAOpModel(core_lep0, LP_SOLVER)
    set_linear_obj!(th_opm, obj_idx, MAX_SENSE)

    # BlobBatch
    bb = BlobBatch(procdir(PROJ, [SIMVER], "batch", (;h=rand(UInt64))))
    # mkpath(bb)
    bb["core_koma"] = Dict[]
    
    for ko in 1:Int(1e8)
        
        # init
        bounds!(th_opm, :, lb0, ub0) # reset bounds
        # th_target_rxn0is = shuffle(target_rxn0is)
        th_target_rxn0is = target_rxn0is
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
                # toko = pop!(th_target_rxn0is)
                toko = rand(th_target_rxn0is)
                l, u = bounds(th_opm, toko)
                bounds!(th_opm, toko, l * DOWNREG_FACTOR, u * DOWNREG_FACTOR)
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
                # _log("ERROR"; err = err_str(e; max_len = 100))
                _break = true;
            end 
            _break && break # for toss
        end # for toss

        # up new koma
        roll_count += 1
        if status != :REVISED
            # push object
            obj = Dict{String, Any}(
                "koset" => koset, 
                "status" => status
            )
            push!(bb["core_koma"], obj)
            # insert hash
            i = searchsortedfirst(koma_hashs, koset_hash)
            insert!(koma_hashs, i, koset_hash)
        end
        
        tot_time = time() - init_time
        effitiency = length(bb["core_koma"]) / roll_count

        # sync blobs
        if length(bb["core_koma"]) > max_nblobs
            
            # up koma_hashs
            _sync_koma_hashs!(koma_hashs)
            
            # up BlobBatch
            bb["meta"]["core_koma.ver"] = ALG_VER # sign
            serialize(bb)

            # empty to save memmory ?
            empty!(bb)
            GC.gc()
            
            # new BlobBatch
            bb = BlobBatch(procdir(PROJ, [SIMVER], "batch", (;h=rand(UInt64))))
            mkpath(bb)
            bb["core_koma"] = Dict[]

            roll_count = 0
        end

        # LOG
        if time() - last_log > log_frec
            lock(PROJ) do
                print("[", getpid(), ".", threadid(), "] ")
                print(" INFO ")
                print("effitiency = ", effitiency)
                print(", koma_hashs_len = ", length(koma_hashs))
                print(", koset_len = ", length(koset))
                print(", effitiency = ", effitiency)
                print(", DOWNREG_FACTOR = ", DOWNREG_FACTOR)
                print(", opt_reltime = ", opt_time / tot_time)
                print(", roll_count = ", roll_count)
                print(", obj_reg_len = ", length(bb["core_koma"]))
                print(", obj_val = ", obj_val)
                print(", batch_size = ", batch_size)
                print(", status = ", status,)
                println()
                println(bb)
            end
            last_log = time()
        end

        if roll_count > (0.1 * max_nblobs) && effitiency < effitiency_th
            break;
        end
    end # for ko

end # @tempcontext