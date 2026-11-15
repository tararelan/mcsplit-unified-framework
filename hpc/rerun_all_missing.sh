#!/bin/bash
# rerun_all_missing.sh - resubmit the missing pairs for EVERY algorithm that has some, for one dataset.   (run on the cluster)
# Usage: bash hpc/rerun_all_missing.sh <mivia|lv|bi> [--swap] [--mem=40G] [--time=..] [--dry-run]
#   Refresh the lists first:  python hpc/find_missing.py
# Algorithms with a reference binary are run with --with-ref, the RL variants without (as in the original runs).
DS=$1
[ -z "$DS" ] && { echo "Usage: bash hpc/rerun_all_missing.sh <mivia|lv|bi> [--swap] [--mem=40G] [--time=..] [--dry-run]"; exit 1; }
OPTS=("${@:2}")
for a in mcsplit down rl ll dal dsb rrsplit symsplit; do
    [ -f "hpc/missing/${a}_${DS}.txt" ] || continue
    bash hpc/resubmit_missing.sh "$a" "$DS" --with-ref "${OPTS[@]}"
done
for a in rl_pair rl_rate rl_down; do
    [ -f "hpc/missing/${a}_${DS}.txt" ] || continue
    bash hpc/resubmit_missing.sh "$a" "$DS" "${OPTS[@]}"
done
