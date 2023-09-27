(julia -t5 --project scripts/BEG2007/CORE_PHASES_KOMA/7_fva.jl -- "SIMVER:ECOLI-CORE-BEG2007-PHASE_I-0.1.0") &
(julia -t5 --project scripts/BEG2007/CORE_PHASES_KOMA/7_shadow_price.jl -- "SIMVER:ECOLI-CORE-BEG2007-PHASE_II-0.1.0") &
(julia -t5 --project scripts/BEG2007/CORE_PHASES_KOMA/2.1_core_xlep0.jl -- "SIMVER:ECOLI-CORE-BEG2007-PHASE_I-0.1.0") &
(julia -t5 --project scripts/BEG2007/CORE_PHASES_KOMA/9_biomass_fba.jl -- "SIMVER:ECOLI-CORE-BEG2007-PHASE_III-0.1.0") &
julia -t5 --project scripts/BEG2007/CORE_PHASES_KOMA/_summary.jl  --  "SIMVER:"

# ============================================================
# ------------------------------------------------------------
# SIMVER: ECOLI-CORE-BEG2007-PHASE_I-0.1.0
# core_biomass_fba.ver: 1066 -- 0.62
# core_strip.ver: 1718 -- 1.0
# batches: 1718 -- 1.0
# core_feasets.ver: 1718 -- 1.0
# core_nut_sp.ver: 1718 -- 1.0
# core_koma.ver: 1718 -- 1.0
# core_fva.ver: 462 -- 0.27
# ------------------------------------------------------------
# SIMVER: ECOLI-CORE-BEG2007-PHASE_III-0.1.0
# core_biomass_fba.ver: 1868 -- 0.7
# core_strip.ver: 2665 -- 1.0
# batches: 2665 -- 1.0
# core_feasets.ver: 2665 -- 1.0
# core_nut_sp.ver: 2655 -- 1.0
# core_koma.ver: 2665 -- 1.0
# core_fva.ver: 1692 -- 0.63
# ------------------------------------------------------------
# SIMVER: ECOLI-CORE-BEG2007-PHASE_II-0.1.0
# core_biomass_fba.ver: 1313 -- 0.63
# core_strip.ver: 2084 -- 1.0
# batches: 2084 -- 1.0
# core_feasets.ver: 2084 -- 1.0
# core_nut_sp.ver: 2078 -- 1.0
# core_koma.ver: 2084 -- 1.0
# core_fva.ver: 1014 -- 0.49


PROCS
2872063 julia -t5 --project scripts/BEG2007/CORE_PHASES_KOMA/6_fva.jl -- SIMVER:ECOLI-CORE-BEG2007-PHASE_II-0.1.0
3030730 julia -t5 --project scripts/BEG2007/CORE_PHASES_KOMA/6_fva.jl -- SIMVER:ECOLI-CORE-BEG2007-PHASE_III-0.1.0
3069266 julia -t5 --project scripts/BEG2007/CORE_PHASES_KOMA/9_biomass_fba.jl -- SIMVER:ECOLI-CORE-BEG2007-PHASE_II-0.1.0
3069405 julia -t5 --project scripts/BEG2007/CORE_PHASES_KOMA/9_biomass_fba.jl -- SIMVER:ECOLI-CORE-BEG2007-PHASE_II-0.1.0
3069552 julia -t5 --project scripts/BEG2007/CORE_PHASES_KOMA/9_biomass_fba.jl -- SIMVER:ECOLI-CORE-BEG2007-PHASE_III-0.1.0
3069558 julia -t5 --project scripts/BEG2007/CORE_PHASES_KOMA/9_biomass_fba.jl -- SIMVER:ECOLI-CORE-BEG2007-PHASE_III-0.1.0
3074721 julia -t5 --project scripts/BEG2007/CORE_PHASES_KOMA/9_biomass_fba.jl -- SIMVER:ECOLI-CORE-BEG2007-PHASE_I-0.1.0
3074752 julia -t5 --project scripts/BEG2007/CORE_PHASES_KOMA/9_biomass_fba.jl -- SIMVER:ECOLI-CORE-BEG2007-PHASE_I-0.1.0