#!/bin/bash
ALGO=$1
A=$2
B=$3

OUTDIR="results/hpc_ref_validation/${ALGO}"
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

# REF_SIZE=""
# REF_NODES=""
# REF_TIME=""

# case "$ALGO" in
#     mcsplit)
#         REF_OUT=$(timeout 1000 ./reference/mcsplit/mcsp -q "$A" "$B" 2>/dev/null)
#         REF_SIZE=$(echo "$REF_OUT" | grep "Solution size" | awk '{print $3}')
#         REF_NODES=$(echo "$REF_OUT" | grep "^Nodes:" | awk '{print $NF}')
#         REF_TIME_MS=$(echo "$REF_OUT" | grep "^CPU time" | awk '{print $NF}')
#         REF_TIME=$(echo "$REF_TIME_MS" | awk '{printf "%.6f", $1/1000}')
#         ;;
#     rl)
#         REF_OUT=$(timeout 1000 ./reference/mcsplit-rl/mcsp+RL min_max -q "$A" "$B" 2>/dev/null)
#         REF_SIZE=$(echo "$REF_OUT" | grep "Solution size" | awk '{print $3}')
#         REF_NODES=$(echo "$REF_OUT" | grep "^Nodes:" | awk '{print $NF}')
#         REF_TIME_MS=$(echo "$REF_OUT" | grep "^CPU time" | awk '{print $NF}')
#         REF_TIME=$(echo "$REF_TIME_MS" | awk '{printf "%.6f", $1/1000}')
#         ;;
#     rllum)
#         REF_OUT=$(timeout 1000 ./reference/mcsplit-rl/mcsp+RL min_max -q "$A" "$B" 2>/dev/null)
#         REF_SIZE=$(echo "$REF_OUT" | grep "Solution size" | awk '{print $3}')
#         REF_NODES=$(echo "$REF_OUT" | grep "^Nodes:" | awk '{print $NF}')
#         REF_TIME_MS=$(echo "$REF_OUT" | grep "^CPU time" | awk '{print $NF}')
#         REF_TIME=$(echo "$REF_TIME_MS" | awk '{printf "%.6f", $1/1000}')
#         ;;
#     dal)
#         REF_OUT=$(timeout 1000 ./reference/mcsplit-dal/mcspDAL min_max -q "$A" "$B" 2>/dev/null)
#         REF_SIZE=$(echo "$REF_OUT" | grep "^#2:" | awk '{print $2}')
#         REF_NODES=$(echo "$REF_OUT" | grep "^#2:" | awk '{print $3}')
#         REF_TIME_MS=$(echo "$REF_OUT" | grep "^#2:" | awk '{print $4}')
#         REF_TIME=$(echo "$REF_TIME_MS" | awk '{printf "%.6f", $1/1000}')
#         ;;
#     ll)
#         REF_OUT=$(timeout 1000 ./reference/mcsplit-ll/mcsp+ll min_max -q "$A" "$B" 2>/dev/null)
#         REF_SIZE=$(echo "$REF_OUT" | grep "Solution size" | awk '{print $3}')
#         REF_NODES=$(echo "$REF_OUT" | grep "^Nodes:" | awk '{print $NF}')
#         REF_TIME_MS=$(echo "$REF_OUT" | grep "^CPU time" | awk '{print $NF}')
#         REF_TIME=$(echo "$REF_TIME_MS" | awk '{printf "%.6f", $1/1000}')
#         ;;
#     dsb)
#         REF_OUT=$(timeout 1000 ./reference/mcsplit-dsb/bin/run.o min_max -q "$A" "$B" 2>/dev/null)
#         REF_SIZE=$(echo "$REF_OUT" | awk -F',' '{print $1}' | tr -d ' ')
#         REF_TIME=$(echo "$REF_OUT" | awk -F',' '{print $4}' | tr -d ' ')
#         REF_NODES=$(echo "$REF_OUT" | awk -F',' '{print $5}' | tr -d ' ')
#         ;;
#     rrsplit)
#         REF_OUT=$(timeout 1000 ./reference/rrsplit/mcsp min_max -q "$A" "$B" 2>/dev/null)
#         REF_SIZE=$(echo "$REF_OUT" | grep "^#2:" | awk '{print $2}')
#         REF_NODES=$(echo "$REF_OUT" | grep "^#2:" | awk '{print $3}')
#         REF_TIME_MS=$(echo "$REF_OUT" | grep "^#2:" | awk '{print $4}')
#         REF_TIME=$(echo "$REF_TIME_MS" | awk '{printf "%.6f", $1/1000}')
#         ;;
#     symsplit)
#         REF_OUT=$(timeout 1000 ./reference/symsplit/bin/run.o min_max -q "$A" "$B" 2>/dev/null)
#         REF_SIZE=$(echo "$REF_OUT" | awk -F',' '{print $1}' | tr -d ' ')
#         REF_TIME=$(echo "$REF_OUT" | awk -F',' '{print $4}' | tr -d ' ')
#         REF_NODES=$(echo "$REF_OUT" | awk -F',' '{print $5}' | tr -d ' ')
#         ;;
# esac

# MATCH="UNKNOWN"
# if [ "$UNIFIED_SIZE" == "$REF_SIZE" ]; then
#     MATCH="PASS"
# elif [ -z "$UNIFIED_SIZE" ] || [ -z "$REF_SIZE" ]; then
#     MATCH="MISSING_DATA"
# else
#     MATCH="FAIL"
# fi

# echo "$(basename "$A"),$(basename "$B"),$ALGO,$UNIFIED_SIZE,$UNIFIED_EDGES,$UNIFIED_NODES,$UNIFIED_TIME,$UNIFIED_ABORTED,$UNIFIED_ROOT_UB,$UNIFIED_NODES_TO_BEST,$UNIFIED_TIME_TO_BEST,$UNIFIED_CUT_BRANCHES,$UNIFIED_BOUND_PRUNED,$UNIFIED_SYM_PRUNED,$UNIFIED_CONFLICTS,$REF_SIZE,$REF_NODES,$REF_TIME,$MATCH" > "$OUTFILE"

echo "$(basename "$A"),$(basename "$B"),$ALGO,$UNIFIED_SIZE,$UNIFIED_EDGES,$UNIFIED_NODES,$UNIFIED_TIME,$UNIFIED_ABORTED,$UNIFIED_ROOT_UB,$UNIFIED_NODES_TO_BEST,$UNIFIED_TIME_TO_BEST,$UNIFIED_CUT_BRANCHES,$UNIFIED_BOUND_PRUNED,$UNIFIED_SYM_PRUNED,$UNIFIED_CONFLICTS" > "$OUTFILE"