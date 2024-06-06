## ------------------------------------------------------------
@time begin
    using Plots
    using Random
    using MetXGEMs
    using MetXBase
    using MetXOptim
    using Statistics
    using MetXNetHub
    using Base.Threads
    using ProgressMeter
    using Combinatorics
    using NutrientLimitedGEMs
end

# TODO: Add Graphs.jl kind of functionality for getting basic stuff
# ------------------------------------------------------------
include("1_setup_sim.jl")
include("1.1_utils.jl")

## ------------------------------------------------------------
@tempcontext ["FVA" => v"0.3.0"] let
    
    # dbs
    glob_db = query(["ROOT", "GLOBALS"])
    xlep_db = query(["ROOT", "XLEP"])

    # koma files
    batch_size = 3 # kos per roll
    objfiles = readdir(procdir(PROJ, [SIMVER]); join = true)
    @threads :static for fn in shuffle(objfiles)
        contains(basename(fn), "obj_reg") || continue
        println("[", getpid(), ".", threadid(), "] ", fn)

        # deserialize
        _, obj_reg = ldat(fn)

        # globals
        LP_SOLVER = glob_db["LP_SOLVER"]
        
        # lep
        elep0 = xlep_db["elep0"][]
        lep0 = lepmodel(elep0)
        M, N = size(lep0)
        
        # run
        do_save = false
        alg_ver = context("FVA")
        for obj in obj_reg

            # check done
            get(obj, "fva_ver", :NONE) == alg_ver && break
            
            koset = obj["koset"]
            feaset = koset[1:(end-batch_size)]
            
            if isempty(feaset) 
                obj["fva_ver"] = alg_ver
                do_save = true
                continue
            end

            # fva
            ko_th = 1e-5
            _with_kos(lep0, feaset) do
                fvalb, fvaub = fva(lep0, LP_SOLVER; verbose = false)
                obj["fva.fvalb"] = fvalb
                obj["fva.fvaub"] = fvaub
                obj["fva.kos"] = findall((abs.(fvalb) .+ abs.(fvaub)) .< ko_th)
                do_save = true
            end 
            obj["fva_ver"] = alg_ver
            
        end # for reg in obj_reg
        
        # serialize
        do_save && sdat(obj_reg, fn)
    end # for fn 
end