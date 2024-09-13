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
    using DataFileNames
    using NutrientLimitedGEMs_2023
end

# ------------------------------------------------------------
include("1_setup_sim.jl")
include("1.1_utils.jl")

## ------------------------------------------------------------
@tempcontext ["NUT_ESSENTIALITY" => v"0.1.0"] let

    # koma files
    batch_size = 3 # kos per roll
    objfiles = readdir(procdir(PROJ, [SIMVER]); join = true)
    @threads :static for fn in objfiles
        contains(basename(fn), "obj_reg") || continue
        @show fn

        # deserialize
        _, obj_reg = ldat(fn)

        # globals
        glob_db = query(["ROOT", "GLOBALS"])
        LP_SOLVER = glob_db["LP_SOLVER"]
        NET_ID = glob_db["NET_ID"]
        NTHREADS = glob_db["NTHREADS"]
        nut_ids = ["EX_glc__D_e", "EX_lac__D_e", "EX_ac_e"]
        
        # xlep
        xlep_db = query(["ROOT", "XLEP"])
        elep0 = xlep_db["elep0"][]
        lep0 = lepmodel(elep0)
        M, N = size(lep0)
        obj_id = extras(lep0, "BIOM")
        obj_idx = colindex(lep0, obj_id)
        alg_ver = context("NUT_ESSENTIALITY")
    
        # run
        for reg in obj_reg

            # check done
            get(reg, "nut_ess_ver", :NONE) == alg_ver && break
            
            koset = reg["koset"]
            feaset = koset[1:(end-batch_size)]
            
            if isempty(feaset) 
                reg["nut_ess_ver"] = alg_ver
                continue
            end

            _with_kos(lep0, feaset) do
                opm = FBAOpModel(lep0, LP_SOLVER)
                set_linear_obj!(opm, obj_idx, MAX_SENSE)
                optimize!(opm)
                obj0 = solution(opm, obj_idx)
                reg["biom"] = Dict{String, Float64}()
                reg["biom"]["feaset"] = obj0
                # println("-"^60)
                # @show obj0

                for ex_id in nut_ids
                    _with_kos(opm, [ex_id]) do 
                        obj = 0.0
                        try; optimize!(opm)
                            obj = solution(opm, obj_idx)
                        catch e; end
                        reg["biom"][ex_id] = obj
                        # @show ex_id
                        # @show 
                    end
                end # for ex_id
            end # _with_kos(lep0, kosets)
            reg["nut_ess_ver"] = alg_ver
        end # for reg in obj_reg
        
        # serialize
        sdat(obj_reg, fn)
    end # for fn 
end

## ------------------------------------------------------------