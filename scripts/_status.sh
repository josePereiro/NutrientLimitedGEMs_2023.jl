while true; do
    clear;

    echo "====================================";
    echo "PROCS";
    pgrep -fla 'ECOLI_CORE_KOMA';

    echo;
    echo "------------------------------------";
    ls data/processed/ECOLI-CORE-0.1.0 -hll | head -n4;

    echo;
    echo "------------------------------------";
    echo;
    tail -n2 data/processed/ECOLI-CORE-0.1.0/koma.log

    echo;

    echo "====================================";
    echo "PROCS";
    pgrep -fla 'iJO1366_KOMA';

    echo;
    echo "------------------------------------";
    ls data/processed/iJO1366-0.1.0 -hll | head -n4;

    echo;
    echo "------------------------------------";
    echo;
    tail -n2 data/processed/iJO1366-0.1.0/koma.log

    echo;
    sleep 5;
    
done