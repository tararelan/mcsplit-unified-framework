#!/bin/bash
# merge_results.sh - concatenate per-instance CSVs into one file with the correct header
# Usage: bash hpc/merge_results.sh <algo> <dataset> [--with-ref] [--swap]
# Use the same flags the run was submitted with (--swap reads the _swap folder).

ALGO=$1
DATASET=$2
WITH_REF=0
SWAP=0
for opt in "${@:3}"; do
    case "$opt" in
        --with-ref) WITH_REF=1 ;;
        --swap)     SWAP=1 ;;
        *) echo "Unknown option: $opt"; exit 1 ;;
    esac
done

if [ -z "$ALGO" ] || [ -z "$DATASET" ]; then
    echo "Usage: bash hpc/merge_results.sh <algo> <dataset> [--with-ref] [--swap]"
    exit 1
fi

if [ "$WITH_REF" = "1" ]; then
    SRCDIR="hpc/results/${ALGO}_${DATASET}_with_ref"
    HEADER="instance_a,instance_b,algo,unified_size,unified_edges,unified_nodes,unified_time,unified_aborted,unified_nodes_to_best,unified_time_to_best,unified_cut_branches,unified_bound_pruned,unified_sym_pruned,unified_bidomain_mean,unified_bidomain_max,ref_size,ref_nodes,ref_time,ref_bidomain_mean,ref_bidomain_max,match"
else
    SRCDIR="hpc/results/${ALGO}_${DATASET}"
    HEADER="instance_a,instance_b,algo,unified_size,unified_edges,unified_nodes,unified_time,unified_aborted,unified_nodes_to_best,unified_time_to_best,unified_cut_branches,unified_bound_pruned,unified_sym_pruned,unified_bidomain_mean,unified_bidomain_max"
fi

# swap runs live in their own folder (see run_instance.sh)
if [ "$SWAP" = "1" ]; then
    SRCDIR="${SRCDIR}_swap"
fi
OUTFILE="${SRCDIR}_all.csv"

if [ ! -d "$SRCDIR" ]; then
    echo "ERROR: $SRCDIR does not exist"
    exit 1
fi

echo "$HEADER" > "$OUTFILE"
cat "$SRCDIR"/*.csv >> "$OUTFILE"

echo "Merged $(ls "$SRCDIR" | wc -l) files into $OUTFILE"