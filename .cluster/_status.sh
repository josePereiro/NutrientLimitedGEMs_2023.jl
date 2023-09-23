while true; do
    clear;

    echo "====================================";
    echo "PROCS";
    pgrep -fla 'CORE_TO_GEM_KOMA';
    pgrep -fla 'ECOLI-CORE-BEG2007';

    # echo;
    # echo "------------------------------------";
    # ls data/processed/ECOLI-CORE-0.1.0 -hll | head -n4;
    # ls data/processed/CORE_TO_GEM-ECOLI-0.1.0 -hll | head -n4;

    echo;
    echo "------------------------------------";
    echo "STORAGE USE"
    du -hs data/processed
    du -hs data/processed/ECOLI-CORE-BEG2007-PHASE_I-0.1.0
    du -hs data/processed/ECOLI-CORE-BEG2007-PHASE_II-0.1.0
    du -hs data/processed/ECOLI-CORE-BEG2007-PHASE_III-0.1.0
    
    echo;
    sleep 5;
    
done