import Pkg

## ---------------------------------------------------
proj = Pkg.project()
proj.name == "NutrientLimitedGEMs" || error("Bad project, ", proj.name)

## ---------------------------------------------------
@info("DEVING")
todev = ["ContextDBs", "ImgTools", "MetX", "MetXBase", "MetXCultureHub", "MetXEP", "MetXNetHub", "MetXOptim", "ProjFlows"]
Pkg.develop(todev)
Pkg.precompile()

## ---------------------------------------------------
import NutrientLimitedGEMs
@info("DONE")