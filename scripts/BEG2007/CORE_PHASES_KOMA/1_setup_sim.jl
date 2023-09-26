using Clp
using ProjFlows
using ContextDBs
using Base.Threads
using NutrientLimitedGEMs

# ------------------------------------------------------------------
# GLOBALS
# ------------------------------------------------------------------

# ------------------------------------------------------------------
#  Project
# ------------------------------------------------------------------

PROJ = Project0(NutrientLimitedGEMs)

# ------------------------------------------------------------------
# ContextDB
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# Utils
function _load_contextdb(id)
    lock(PROJ) do
        name = string("db.", id, ".jls")
        fn, DB = withprocdat(PROJ, :get!, name) do 
            ContextDB()
        end
        contextdb!(DB)
        println(fn)
    end
end

function _save_contextdb(id)
    lock(PROJ) do
        name = string("db.", id, ".jls")
        fn, _ = withprocdat(PROJ, :set!, name) do 
            contextdb()
        end
        println(fn)
    end
end

# ------------------------------------------------------------------
# ARGS
if isinteractive()
    # ARGS DEV
    SIMVER = "ECOLI-CORE-BEG2007-PHASE_I-0.1.0"
else
    SIMVER = parseARGS("SIMVER:")
end

# ------------------------------------------------------------------
_load_contextdb(SIMVER)
cacherefs_dir!(cachedir(PROJ, SIMVER))
@context! "ROOT" SIMVER
@commit! ["META"] "DESC" => """
    # Reproducing the enviroment of phase 2 (only glc) from Beg 2007
    # 1. Using the core model (reduced)
    #     1.1. I'll impose fractional KO (flux reduction) batches to the network and evaluate:
    #         i) feasibility
    #         ii) shadow price for several nutrients
    #         iii) etc
    # 2. Using the GEM (iJO1366)
    #     2.1 Use the core samples as template to produce GEM samples, (a projection) and evaluate:
    #         i) feasibility 
    #         ii) shadow price for several nutrients
    #         iii) etc
    # 3. Compare the two results with data from Q. K. Beg, 2007 testing at least the separation glucose from the rest
"""

# ------------------------------------------------------------------
@tempcontext ["GLOBALS"] begin
    @stage! CORE_NET_ID = "ecoli_core"
    @stage! GEM_NET_ID = "iJO1366"
    @stage! LP_SOLVER = Clp.Optimizer
    @stage! NTHREADS = max(nthreads() - 2, 1)
end

# ------------------------------------------------------------------
return nothing