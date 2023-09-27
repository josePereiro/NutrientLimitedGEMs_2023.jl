using Clp
using ProjFlows
using ContextDBs
using Base.Threads
using NutrientLimitedGEMs

# ------------------------------------------------------------------
#  Project
# ------------------------------------------------------------------

PROJ = Project0(NutrientLimitedGEMs)

# ------------------------------------------------------------------
include("1.1_utils.jl")

# ------------------------------------------------------------------
# ContexDB
# ------------------------------------------------------------------
SIMVER = "ECOLI-CORE-BEG2007-MICROARRAY-0.1.0"

# ------------------------------------------------------------------
_load_contextdb(SIMVER)
cacherefs_dir!(cachedir(PROJ, SIMVER))
@context! "ROOT" SIMVER
@commit! ["META"] "DESC" => """
    Data from Beg 2007 Supp1 and Supp3 is collected and matched against 
    downregulation ensembles from ecoli core. 
    We can evaluate ensembles by its similarity with regulatory data.
"""

# GLOBALS
@tempcontext ["GLOBALS"] begin
    @stage! CORE_NET_ID = "ecoli_core"
    @stage! GEM_NET_ID = "iJO1366"
    @stage! LP_SOLVER = Clp.Optimizer
    @stage! NTHREADS = max(nthreads() - 2, 1)
    @stage! DOWNREG_FACTOR = 0.3
    @stage! KO_OBJ_VAL_TH = 0.01
    @stage! DOWNREG_BATCH_SIZE = 3
    @stage! CORE_KOMA_PHASES = [
        "ECOLI-CORE-BEG2007-PHASE_I-0.1.0",
        "ECOLI-CORE-BEG2007-PHASE_II-0.1.0",
        "ECOLI-CORE-BEG2007-PHASE_III-0.1.0",
    ] # TOSYNC
end

# ------------------------------------------------------------------
return nothing