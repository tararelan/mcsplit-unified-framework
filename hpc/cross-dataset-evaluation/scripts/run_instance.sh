#!/bin/bash
ALGO=$1
A=$2
B=$3
DATASET=$4  # bi or lv

OUTDIR="hpc/cross-dataset-evaluation/results/${ALGO}_${DATASET}"
mkdir -p "$OUTDIR"

OUTFILE="${OUTDIR}/$(basename "$A").$(basename "$B").csv"

if [ -f "$OUTFILE" ]; then
    exit 0
fi

FORMAT_FLAG=""
if [ "$DATASET" = "bi" ]; then
    FORMAT_FLAG="-l -i"
elif [ "$DATASET" = "lv" ]; then
    FORMAT_FLAG="-l"
fi

# NOTE: fixed a pre-existing off-by-one bug here - this used to read fields
# $7-$11 for nodes_to_best..sym_pruned, skipping $6 entirely, which meant
# every one of those five columns was actually reading the NEXT field over
# (e.g. "nodes_to_best" was really time_to_best). Corrected to $6-$10, and
# added the bidomain mean/max fields (11-12) the binary now also prints.
UNIFIED_RESULT=$(./solvers/bin/mcsp -A "$ALGO" -t 1000 $FORMAT_FLAG -q "$A" "$B" 2>/dev/null)
UNIFIED_SIZE=$(echo "$UNIFIED_RESULT" | awk '{print $1}')
UNIFIED_EDGES=$(echo "$UNIFIED_RESULT" | awk '{print $2}')
UNIFIED_NODES=$(echo "$UNIFIED_RESULT" | awk '{print $3}')
UNIFIED_TIME=$(echo "$UNIFIED_RESULT" | awk '{print $4}')
UNIFIED_ABORTED=$(echo "$UNIFIED_RESULT" | awk '{print $5}')
UNIFIED_NODES_TO_BEST=$(echo "$UNIFIED_RESULT" | awk '{print $6}')
UNIFIED_TIME_TO_BEST=$(echo "$UNIFIED_RESULT" | awk '{print $7}')
UNIFIED_CUT_BRANCHES=$(echo "$UNIFIED_RESULT" | awk '{print $8}')
UNIFIED_BOUND_PRUNED=$(echo "$UNIFIED_RESULT" | awk '{print $9}')
UNIFIED_SYM_PRUNED=$(echo "$UNIFIED_RESULT" | awk '{print $10}')
UNIFIED_BIDOMAIN_MEAN=$(echo "$UNIFIED_RESULT" | awk '{print $11}')
UNIFIED_BIDOMAIN_MAX=$(echo "$UNIFIED_RESULT" | awk '{print $12}')

echo "$(basename "$A"),$(basename "$B"),$ALGO,$UNIFIED_SIZE,$UNIFIED_EDGES,$UNIFIED_NODES,$UNIFIED_TIME,$UNIFIED_ABORTED,$UNIFIED_NODES_TO_BEST,$UNIFIED_TIME_TO_BEST,$UNIFIED_CUT_BRANCHES,$UNIFIED_BOUND_PRUNED,$UNIFIED_SYM_PRUNED,$UNIFIED_BIDOMAIN_MEAN,$UNIFIED_BIDOMAIN_MAX" > "$OUTFILE"