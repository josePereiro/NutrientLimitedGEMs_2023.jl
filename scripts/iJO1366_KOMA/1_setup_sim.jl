using Gurobi
using ProjFlows
using ContextDBs
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
    name = string("db.", id, ".jls")
    fn, DB = withprocdat(PROJ, :get!, name) do 
        ContextDB()
    end
    contextdb!(DB)
    println(fn)
end

function _save_contextdb(id)
    name = string("db.", id, ".jls")
    fn, _ = withprocdat(PROJ, :set!, name) do 
        contextdb()
    end
    println(fn)
end

# ------------------------------------------------------------------
# Top context

SIMVER = "iJO1366-0.1.0"

_load_contextdb(SIMVER)
cacherefs_dir!(cachedir(PROJ))
@context! "ROOT" SIMVER
@commit! ["META"] "DESC" => """
    1. I'll impose KO batches to the network and evaluate:
    i) feasibility 
    ii) shadow price for several nutrients
"""

# ------------------------------------------------------------------
@tempcontext ["GLOBALS"] begin
    @stage! NET_ID = "iJO1366"
    @stage! LP_SOLVER = Gurobi.Optimizer
end


# ------------------------------------------------------------------
return nothing