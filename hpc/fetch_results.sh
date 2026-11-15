#!/bin/bash
# fetch_results.sh - copy the merged result files of one dataset from the cluster into hpc/publication.   (run on YOUR computer, e.g. in WSL, from the repo root)
# Usage: bash hpc/fetch_results.sh <mivia|lv|bi>
# The current published files are first copied to hpc/publication/backup_before_rerun/ (existing backups are not overwritten).
DS=$1
HOST=tr77@hypatia.st-andrews.ac.uk
REMOTE=work/tr77-cs5099/hpc/results
[ -z "$DS" ] && { echo "Usage: bash hpc/fetch_results.sh <mivia|lv|bi>"; exit 1; }
mkdir -p hpc/publication/backup_before_rerun
for f in hpc/publication/*_${DS}_*all.csv; do cp -n "$f" hpc/publication/backup_before_rerun/; done
for a in mcsplit down rl ll dal dsb rrsplit symsplit; do
    scp "$HOST:$REMOTE/${a}_${DS}_with_ref_all.csv" "hpc/publication/${a}_${DS}_with_ref_all.csv"
done
for a in rl_pair rl_rate rl_down; do
    scp "$HOST:$REMOTE/${a}_${DS}_all.csv" "hpc/publication/${a}_${DS}_all.csv"
done
