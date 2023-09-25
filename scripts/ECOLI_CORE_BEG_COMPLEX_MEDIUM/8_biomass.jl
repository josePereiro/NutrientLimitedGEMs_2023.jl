## ------------------------------------------------------------
@time begin
    using Random
    using MetXGEMs
    using MetXNetHub
    using MetXBase
    using MetXOptim
    using Base.Threads
    using NutrientLimitedGEMs
end

# TODO: Add Graphs.jl kind of functionality for getting basic stuff
# ------------------------------------------------------------
include("1_setup_sim.jl")
include("1.1_utils.jl")

## ------------------------------------------------------------
@tempcontext ["CORE_BIOMASS" => v"0.1.0"] let
    
    # dbs
    glob_db = query(["ROOT", "GLOBALS"])
    xlep_db = query(["ROOT", "CORE_XLEP"])

    # koma files
    downreg_factor = 0.3 # TOSYNC
    objfiles = readdir(procdir(PROJ, [SIMVER]); join = true)
    # @threads 
    for fn in shuffle(objfiles)
        contains(basename(fn), "obj_reg") || continue

        # deserialize
        @time obj_reg = try_ldat(fn)
        isempty(obj_reg) && continue

        # globals
        LP_SOLVER = glob_db["LP_SOLVER"]
        
        # lep
        core_elep0 = xlep_db["core_elep0"][]
        core_lep0 = lepmodel(core_elep0)
        biom_id = extras(core_lep0, "BIOM")
        core_elep0 = nothing

        # opm
        opm = FBAOpModel(core_lep0, LP_SOLVER)
        
        # run
        do_save = false
        ALG_VER = context("CORE_BIOMASS")
        info_frec = 100
        gc_frec = 10
        for (obji, obj) in enumerate(obj_reg)
            
            # check done
            get(obj, "core_biomass.ver", :NONE) == ALG_VER && continue
            haskey(obj, "core_strip.koset") || continue
            haskey(obj, "core_feasets") || continue
            haskey(obj, "core_nut_sp.ver") || continue
            do_save = true

            # info
            show_flag = obji == 1 || obji == lastindex(obj_reg) || iszero(rem(obji, info_frec)) 
            show_flag && println("[", getpid(), ".", threadid(), "] ", 
                "obji ", obji, "\\", length(obj_reg), " ",
                basename(fn)
            )
            
            # compute
            koset = obj["core_strip.koset"]
            feasets = obj["core_feasets"]
            for (li, feaobj) in feasets
                feaset = koset[1:li]
                _with_downreg(opm, feaset, downreg_factor) do
                    optimize!(opm)
                    feaobj["core_biomass.biom"] = solution(opm, biom_id)
                end  # _with_downreg
            end # for feasets

            # ALG_VER
            obj["core_biomass.ver"] = ALG_VER

            # GC (TEST)
            gc_flag = iszero(rem(obji, gc_frec))
            gc_flag && GC.gc()

        end # for reg in obj_reg
        
        # serialize
        do_save && sdat(obj_reg, fn)

        # TEST
        do_save && exit()

    end # for fn 
end