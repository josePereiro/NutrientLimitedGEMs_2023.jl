while true; do
    cd /home/pereiro/dev/NutrientLimitedGEMs
    pwd
    julia --project -t30 scripts/CORE_TO_GEM_KOMA/4_strip.jl || true
    julia --project -t30 scripts/CORE_TO_GEM_KOMA/5_fva.jl || true
    sleep 10
done

(julia --project -t5 scripts/CORE_TO_GEM_KOMA/4_strip.jl) &