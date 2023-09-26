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
@tempcontext ["GLOBALS"] begin
    @stage! SIMVER = "ECOLI-CORE-BEG2007-MICROARRAY-0.1.0"
    @stage! CORE_NET_ID = "ecoli_core"
    @stage! GEM_NET_ID = "iJO1366"
    @stage! LP_SOLVER = Clp.Optimizer
    @stage! NTHREADS = max(nthreads() - 2, 1)
    @stage! DOWNREG_FACTOR = 0.3
    @stage! KO_OBJ_VAL_TH = 0.01
    @stage! DOWNREG_BATCH_SIZE = 3
end

# ------------------------------------------------------------------
_load_contextdb(SIMVER)
cacherefs_dir!(cachedir(PROJ, SIMVER))
@context! "ROOT" SIMVER
@commit! ["META"] "DESC" => """
    Data from Beg 2007 Supp1 and Supp3 is collected and matched against 
    downregulation ensembles from ecoli core. 
    We can evaluate ensembles by its similarity with regulatory data.
"""

# ------------------------------------------------------------------
return nothing