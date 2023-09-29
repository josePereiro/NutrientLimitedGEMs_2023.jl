(julia -t2 --project scripts/BEG2007/PHASES_KOMA/7_core_fva.jl -- "SIMVER:ECOLI-CORE-BEG2007-PHASE_I-0.1.0") &
(julia -t5 --project scripts/BEG2007/PHASES_KOMA/7_shadow_price.jl -- "SIMVER:ECOLI-CORE-BEG2007-PHASE_II-0.1.0") &
(julia -t5 --project scripts/BEG2007/PHASES_KOMA/2.1_core_xlep0.jl -- "SIMVER:ECOLI-CORE-BEG2007-PHASE_I-0.1.0") &
(julia -t5 --project scripts/BEG2007/PHASES_KOMA/3.3_rxn_map.jl -- "SIMVER:ECOLI-CORE-BEG2007-PHASE_I-0.1.0") &
(julia -t5 --project scripts/BEG2007/PHASES_KOMA/3.2_gem_xlep0.jl -- "SIMVER:ECOLI-CORE-BEG2007-PHASE_I-0.1.0") &
(julia -t5 --project scripts/BEG2007/PHASES_KOMA/9_gem_biomass_fba.jl -- "SIMVER:ECOLI-CORE-BEG2007-PHASE_I-0.1.0") &
(julia -t5 --project scripts/BEG2007/PHASES_KOMA/9_biomass_fba.jl -- "SIMVER:ECOLI-CORE-BEG2007-PHASE_III-0.1.0") &
(julia -t5 --project scripts/BEG2007/PHASES_KOMA/10_core_ep.jl -- "SIMVER:ECOLI-CORE-BEG2007-PHASE_I-0.1.0") &
julia -t5 --project scripts/BEG2007/PHASES_KOMA/_summary.jl  --  "SIMVER:"
(julia -t5 --project scripts/BEG2007/PHASES_KOMA/6.2_plots.jl  --  "SIMVER:") &

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
# ============================================================
# ------------------------------------------------------------
# SIMVER: ECOLI-CORE-BEG2007-PHASE_I-0.1.0
# core_fva.ver: 670 -- 0.39
# core_strip.ver: 1718 -- 1.0
# gem_biomass_fba.ver: 1140 -- 0.66
# batches: 1718 -- 1.0
# core_feasets.ver: 1718 -- 1.0
# core_nut_sp.ver: 1718 -- 1.0
# core_koma.ver: 1718 -- 1.0
# core_biomass_fba.ver: 1718 -- 1.0
# ------------------------------------------------------------
# SIMVER: ECOLI-CORE-BEG2007-PHASE_III-0.1.0
# core_biomass_fba.ver: 2663 -- 1.0
# core_strip.ver: 2665 -- 1.0
# batches: 2665 -- 1.0
# core_feasets.ver: 2665 -- 1.0
# core_nut_sp.ver: 2655 -- 1.0
# core_koma.ver: 2665 -- 1.0
# core_fva.ver: 1717 -- 0.64
# ------------------------------------------------------------
# SIMVER: ECOLI-CORE-BEG2007-PHASE_II-0.1.0
# core_biomass_fba.ver: 2083 -- 1.0
# core_strip.ver: 2084 -- 1.0
# batches: 2084 -- 1.0
# core_feasets.ver: 2084 -- 1.0
# core_nut_sp.ver: 2078 -- 1.0
# core_koma.ver: 2084 -- 1.0
# core_fva.ver: 1024 -- 0.49




