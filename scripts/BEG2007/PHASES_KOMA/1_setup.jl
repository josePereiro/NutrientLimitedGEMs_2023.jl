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

# ARGS
if isinteractive()
    # ARGS DEV
    SIMVER = "ECOLI-CORE-BEG2007-PHASE_I-0.1.0"
else
    SIMVER = parseARGS("SIMVER:")
end

_load_contextdb(SIMVER)
cacherefs_dir!(cachedir(PROJ, SIMVER))
@context! "ROOT" SIMVER
@commit! ["META"] "DESC" => """
    1. Contextualize the ecoli core network to model the 
    phasess of mixed carbon cultures described at Beg et al. 2007 (https://doi.org/10.1073/pnas.0609845104.).
    2. Compute ensamble of koma (downregulation sets which kill the network)
    3. Re sample the kome to get feasible ensemble (downregulation sets with non zero growth)
    4. Compute properties of the ensembles
        - max biomass fba solution
        - fva bounds
        - nutrinet shadow price
        - EP entropy and free energy
        - 
"""

@tempcontext ["GLOBALS"] begin
    @stage! CORE_NET_ID = "ecoli_core"
    @stage! GEM_NET_ID = "iJO1366"
    @stage! LP_SOLVER = Clp.Optimizer
    @stage! NTHREADS = max(nthreads() - 2, 1)
    @stage! DOWNREG_FACTOR = 0.3
    @stage! KO_OBJ_VAL_TH = 0.01
    @stage! DOWNREG_BATCH_SIZE = 3
end

# ------------------------------------------------------------------
return nothing