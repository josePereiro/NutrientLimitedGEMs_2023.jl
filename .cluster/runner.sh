(julia -t5 --project scripts/BEG2007/CORE_PHASES_KOMA/6_fva.jl -- "SIMVER:ECOLI-CORE-BEG2007-PHASE_I-0.1.0") &
(julia -t5 --project scripts/BEG2007/CORE_PHASES_KOMA/7_shadow_price.jl -- "SIMVER:ECOLI-CORE-BEG2007-PHASE_II-0.1.0") &
(julia -t5 --project scripts/BEG2007/CORE_PHASES_KOMA/2.1_core_xlep0.jl -- "SIMVER:ECOLI-CORE-BEG2007-PHASE_I-0.1.0") &
julia -t5 --project scripts/BEG2007/CORE_PHASES_KOMA/_summary.jl  --  "SIMVER:"

