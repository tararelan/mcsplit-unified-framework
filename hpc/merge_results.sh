#!/bin/bash
# merge_results.sh - concatenate per-instance CSVs into one file with the correct header
# Usage: bash hpc/merge_results.sh <algo> <dataset> [--with-ref]

ALGO=$1
DATASET=$2
WITH_REF=0
if [ "$3" = "--with-ref" ]; then
    WITH_REF=1
fi

if [ -z "$ALGO" ] || [ -z "$DATASET" ]; then
    echo "Usage: bash hpc/merge_results.sh <algo> <dataset> [--with-ref]"
    exit 1
fi

if [ "$WITH_REF" = "1" ]; then
    SRCDIR="hpc/results/${ALGO}_${DATASET}_with_ref"
    OUTFILE="hpc/results/${ALGO}_${DATASET}_with_ref_all.csv"
    HEADER="instance_a,instance_b,algo,unified_size,unified_edges,unified_nodes,unified_time,unified_aborted,unified_nodes_to_best,unified_time_to_best,unified_cut_branches,unified_bound_pruned,unified_sym_pruned,ref_size,ref_nodes,ref_time,match"
else
    SRCDIR="hpc/results/${ALGO}_${DATASET}"
    OUTFILE="hpc/results/${ALGO}_${DATASET}_all.csv"
    HEADER="instance_a,instance_b,algo,unified_size,unified_edges,unified_nodes,unified_time,unified_aborted,unified_nodes_to_best,unified_time_to_best,unified_cut_branches,unified_bound_pruned,unified_sym_pruned"
fi

if [ ! -d "$SRCDIR" ]; then
    echo "ERROR: $SRCDIR does not exist"
    exit 1
fi

echo "$HEADER" > "$OUTFILE"
cat "$SRCDIR"/*.csv >> "$OUTFILE"

echo "Merged $(ls "$SRCDIR" | wc -l) files into $OUTFILE"