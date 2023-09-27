(julia -t5 --project scripts/BEG2007/PHASES_KOMA/7_core_fva.jl -- "SIMVER:ECOLI-CORE-BEG2007-PHASE_I-0.1.0") &
(julia -t5 --project scripts/BEG2007/PHASES_KOMA/7_shadow_price.jl -- "SIMVER:ECOLI-CORE-BEG2007-PHASE_II-0.1.0") &
(julia -t5 --project scripts/BEG2007/PHASES_KOMA/2.1_core_xlep0.jl -- "SIMVER:ECOLI-CORE-BEG2007-PHASE_I-0.1.0") &
(julia -t5 --project scripts/BEG2007/PHASES_KOMA/9_biomass_fba.jl -- "SIMVER:ECOLI-CORE-BEG2007-PHASE_III-0.1.0") &
julia -t5 --project scripts/BEG2007/PHASES_KOMA/_summary.jl  --  "SIMVER:"
(julia -t5 --project scripts/BEG2007/PHASES_KOMA/6.2_plots.jl  --  "SIMVER:") &

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
