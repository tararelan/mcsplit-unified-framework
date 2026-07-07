#!/bin/bash
ALGO=$1
A=$2
B=$3

OUTDIR="ablation/results/${ALGO}"
mkdir -p "$OUTDIR"

OUTFILE="${OUTDIR}/$(basename "$A").$(basename "$B").csv"

if [ -f "$OUTFILE" ]; then
    exit 0
fi

UNIFIED_RESULT=$(./solvers/bin/mcsp -A "$ALGO" -t 1000 -q "$A" "$B" 2>/dev/null)
UNIFIED_SIZE=$(echo "$UNIFIED_RESULT" | awk '{print $1}')
UNIFIED_EDGES=$(echo "$UNIFIED_RESULT" | awk '{print $2}')
UNIFIED_NODES=$(echo "$UNIFIED_RESULT" | awk '{print $3}')
UNIFIED_TIME=$(echo "$UNIFIED_RESULT" | awk '{print $4}')
UNIFIED_ABORTED=$(echo "$UNIFIED_RESULT" | awk '{print $5}')
UNIFIED_ROOT_UB=$(echo "$UNIFIED_RESULT" | awk '{print $6}')
UNIFIED_NODES_TO_BEST=$(echo "$UNIFIED_RESULT" | awk '{print $7}')
UNIFIED_TIME_TO_BEST=$(echo "$UNIFIED_RESULT" | awk '{print $8}')
UNIFIED_CUT_BRANCHES=$(echo "$UNIFIED_RESULT" | awk '{print $9}')
UNIFIED_BOUND_PRUNED=$(echo "$UNIFIED_RESULT" | awk '{print $10}')
UNIFIED_SYM_PRUNED=$(echo "$UNIFIED_RESULT" | awk '{print $11}')
UNIFIED_CONFLICTS=$(echo "$UNIFIED_RESULT" | awk '{print $12}')


echo "$(basename "$A"),$(basename "$B"),$ALGO,$UNIFIED_SIZE,$UNIFIED_EDGES,$UNIFIED_NODES,$UNIFIED_TIME,$UNIFIED_ABORTED,$UNIFIED_ROOT_UB,$UNIFIED_NODES_TO_BEST,$UNIFIED_TIME_TO_BEST,$UNIFIED_CUT_BRANCHES,$UNIFIED_BOUND_PRUNED,$UNIFIED_SYM_PRUNED,$UNIFIED_CONFLICTS" > "$OUTFILE"