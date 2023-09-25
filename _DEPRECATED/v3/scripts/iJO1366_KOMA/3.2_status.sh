while true; do
    clear;
    echo "------------------------------------";
    echo "PROCS";
    pgrep -fla 'iJO1366_KOMA';

    echo;
    echo "------------------------------------";
    ls data/processed/ECOLI-CORE-0.1.0 -hll | head;

    echo;
    echo "------------------------------------";
    echo;
    tail data/processed/ECOLI-CORE-0.1.0/koma.log

    echo;
    sleep 5;
done