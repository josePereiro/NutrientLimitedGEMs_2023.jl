(julia -t5 --project scripts/BEG2007/PHASES_KOMA/3.1_core_xlep0.jl -- "SIMVER:ECOLI-CORE-BEG2007-PHASE_0") &
(julia -t7 --project scripts/BEG2007/PHASES_KOMA/3.2_gem_xlep0.jl -- "SIMVER:ECOLI-CORE-BEG2007-PHASE_0") &
(julia -t5 --project scripts/BEG2007/PHASES_KOMA/3.3_rxn_map.jl -- "SIMVER:ECOLI-CORE-BEG2007-PHASE_0") &
(julia -t2 --project scripts/BEG2007/PHASES_KOMA/4_core_koma.jl -- "SIMVER:ECOLI-CORE-BEG2007-PHASE_1") &
(julia -t2 --project scripts/BEG2007/PHASES_KOMA/5_koma_strip.jl -- "SIMVER:ECOLI-CORE-BEG2007-PHASE_1") &
(julia -t5 --project scripts/BEG2007/PHASES_KOMA/6_core_feasets.jl -- "SIMVER:ECOLI-CORE-BEG2007-PHASE_1") &
(julia -t5 --project scripts/BEG2007/PHASES_KOMA/6.1_enumerate.jl -- "SIMVER:ECOLI-CORE-BEG2007-PHASE_1") &
(julia -t2 --project scripts/BEG2007/PHASES_KOMA/7_core_fva.jl -- "SIMVER:ECOLI-CORE-BEG2007-PHASE_1") &
(julia -t5 --project scripts/BEG2007/PHASES_KOMA/7_shadow_price.jl -- "SIMVER:ECOLI-CORE-BEG2007-PHASE_2") &
(julia -t5 --project scripts/BEG2007/PHASES_KOMA/9_gem_biomass_fba.jl -- "SIMVER:ECOLI-CORE-BEG2007-PHASE_1") &
(julia -t5 --project scripts/BEG2007/PHASES_KOMA/9_core_biomass_fba.jl -- "SIMVER:ECOLI-CORE-BEG2007-PHASE_1") &
(julia -t5 --project scripts/BEG2007/PHASES_KOMA/10_core_ep.jl -- "SIMVER:ECOLI-CORE-BEG2007-PHASE_1") &
(julia -t10 --project scripts/BEG2007/PHASES_KOMA/12.1_ensembles_ph1_v1.jl -- "SIMVER:") &

julia -t3 --project scripts/BEG2007/PHASES_KOMA/_summary.jl  --  "SIMVER:"
(julia -t5 --project scripts/BEG2007/PHASES_KOMA/6.2_plots.jl  --  "SIMVER:") &

# Ensembles
(julia -t5 --project scripts/BEG2007/PHASES_KOMA/12.1_ensem_ph1_zU.jl -- "SIMVER:" "ENS-SIZE:5000") &
(julia -t5 --project scripts/BEG2007/PHASES_KOMA/12.1_ensem_zU.jl -- "SIMVER:ECOLI-CORE-BEG2007-PHASE_1" "ENS-SIZE:5000" "BIOMAS-DIST:Uniform") &
(julia -t5 --project scripts/BEG2007/PHASES_KOMA/12.1_ensem_zU.jl -- "SIMVER:ECOLI-CORE-BEG2007-PHASE_2" "ENS-SIZE:5000" "BIOMAS-DIST:Uniform") &
(julia -t5 --project scripts/BEG2007/PHASES_KOMA/12.1_ensem_zU.jl -- "SIMVER:ECOLI-CORE-BEG2007-PHASE_3" "ENS-SIZE:5000" "BIOMAS-DIST:Uniform") &

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
============================================================
------------------------------------------------------------
SIMVER: ECOLI-CORE-BEG2007-PHASE_0
core_biomass_fba.ver: 1718 -- 1.0
core_strip.ver: 1718 -- 1.0
gem_biomass_fba.ver: 1140 -- 0.66
batches: 1718 -- 1.0
core_feasets.ver: 1718 -- 1.0
core_ep.ver: 41 -- 0.024
core_nut_sp.ver: 1718 -- 1.0
core_koma.ver: 1718 -- 1.0
core_fva.ver: 676 -- 0.39
------------------------------------------------------------
SIMVER: ECOLI-CORE-BEG2007-PHASE_1
core_strip.ver: 994 -- 0.73
batches: 1356 -- 1.0
gem_biomass_fba.ver: 88 -- 0.065
core_feasets.ver: 994 -- 0.73
core_koma.ver: 1356 -- 1.0
core_biomass_fba.ver: 376 -- 0.28
------------------------------------------------------------
SIMVER: ECOLI-CORE-BEG2007-PHASE_3
core_biomass_fba.ver: 2663 -- 1.0
core_strip.ver: 2665 -- 1.0
batches: 2665 -- 1.0
core_feasets.ver: 2665 -- 1.0
core_nut_sp.ver: 2655 -- 1.0
core_koma.ver: 2665 -- 1.0
core_fva.ver: 1717 -- 0.64
------------------------------------------------------------
SIMVER: ECOLI-CORE-BEG2007-PHASE_2
core_biomass_fba.ver: 2083 -- 1.0
core_strip.ver: 2084 -- 1.0
batches: 2084 -- 1.0
core_feasets.ver: 2084 -- 1.0
core_nut_sp.ver: 2078 -- 1.0
core_koma.ver: 2084 -- 1.0
core_fva.ver: 1024 -- 0.49
