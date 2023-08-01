import Pkg

## ---------------------------------------------------
proj = Pkg.project()
proj.name == "NutrientLimitedGEMs" || error("Bad project, ", proj.name)

## ---------------------------------------------------
todev = ["ContextDBs", "ImgTools", "MetX", "MetXBase", "MetXCultureHub", "MetXEP", "MetXNetHub", "MetXOptim", "ProjFlows"]
Pkg.develop(todev)
Pkg.precompile()