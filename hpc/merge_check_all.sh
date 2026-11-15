#!/bin/bash
# merge_check_all.sh - merge the result folders of every algorithm for one dataset and report rows / blank rows.   (run on the cluster)
# Usage: bash hpc/merge_check_all.sh <mivia|lv|bi> [--swap]
# Run it only when the jobs have finished. Expected rows: bi 9180, lv 6105, mivia 8140; blank must be 0.
DS=$1
[ -z "$DS" ] && { echo "Usage: bash hpc/merge_check_all.sh <mivia|lv|bi>"; exit 1; }
SW=""; [ "$2" = "--swap" ] && SW="_swap"
running=$(squeue -u "$USER" -h -r | wc -l)
[ "$running" -gt 0 ] && echo "WARNING: $running of your jobs are still queued/running - merge again when they finish"
report() {   # $1 = merged csv
    echo "$(basename "$1"): rows=$(($(wc -l < "$1")-1)) blank=$(awk -F, 'NR>1 && $4==""' "$1" | wc -l)"
}
for a in mcsplit down rl ll dal dsb rrsplit symsplit; do
    [ -d "hpc/results/${a}_${DS}_with_ref${SW}" ] || continue
    bash hpc/merge_results.sh "$a" "$DS" --with-ref ${SW:+--swap} >/dev/null && report "hpc/results/${a}_${DS}_with_ref${SW}_all.csv"
done
for a in rl_pair rl_rate rl_down; do
    [ -d "hpc/results/${a}_${DS}${SW}" ] || continue
    bash hpc/merge_results.sh "$a" "$DS" ${SW:+--swap} >/dev/null && report "hpc/results/${a}_${DS}${SW}_all.csv"
done
