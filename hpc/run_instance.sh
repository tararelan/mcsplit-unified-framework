#!/bin/bash
ALGO=$1
A=$2
B=$3
DATASET=${4:-"mivia"}  # mivia, bi, or lv (default: mivia)
WITH_REF=${5:-0}       # 1 = also run reference binary

# Output directory
if [ "$WITH_REF" = "1" ]; then
    OUTDIR="hpc/results/${ALGO}_${DATASET}_with_ref"
else
    OUTDIR="hpc/results/${ALGO}_${DATASET}"
fi
mkdir -p "$OUTDIR"

OUTFILE="${OUTDIR}/$(basename "$A").$(basename "$B").csv"

if [ -f "$OUTFILE" ]; then
    exit 0
fi

# Format flag: LAD for BI/LV, binary for MIVIA
FORMAT_FLAG=""
if [ "$DATASET" = "bi" ]; then
    FORMAT_FLAG="-l -i"
elif [ "$DATASET" = "lv" ]; then
    FORMAT_FLAG="-l"
fi

# Run unified solver
# Bidomain-count summary (mean, max across the whole search) rides along as
# fields 11-12 of this same summary line, rather than the full per-node
# sequence or a separate dump file - kept small regardless of node count.
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

# Run reference binary (optional)
REF_SIZE=""
REF_NODES=""
REF_TIME=""
REF_BIDOMAIN_MEAN=""
REF_BIDOMAIN_MAX=""

if [ "$WITH_REF" = "1" ]; then
    case "$ALGO" in
        mcsplit)
            # The mcsplit reference binary uses -r (not -i) to mark a graph as
            # directed, unlike every other reference binary and the unified
            # framework - swap it here so directed (bi) runs actually get
            # treated as directed instead of failing on an unrecognised flag.
            MCSPLIT_FORMAT_FLAG=$(echo "$FORMAT_FLAG" | sed 's/-i/-r/')
            REF_OUT=$(timeout 1000 ./reference/mcsplit/mcsp $MCSPLIT_FORMAT_FLAG -q "$A" "$B" 2>/dev/null)
            REF_SIZE=$(echo "$REF_OUT" | grep "Solution size" | awk '{print $3}')
            REF_NODES=$(echo "$REF_OUT" | grep "^Nodes:" | awk '{print $NF}')
            REF_TIME_MS=$(echo "$REF_OUT" | grep "^CPU time" | awk '{print $NF}')
            REF_TIME=$(echo "$REF_TIME_MS" | awk '{printf "%.6f", $1/1000}')
            REF_BIDOMAIN_RAW=$(echo "$REF_OUT" | grep "^BIDOMAINS:" | sed 's/^BIDOMAINS://')
            REF_BIDOMAIN_MEAN=$(echo "$REF_BIDOMAIN_RAW" | cut -d',' -f1)
            REF_BIDOMAIN_MAX=$(echo "$REF_BIDOMAIN_RAW" | cut -d',' -f2)
            ;;
        down)
            # McSplit-down has no separate reference: it is reference/mcsplit/mcsp
            # run with --big-first (-b). Same -i -> -r swap as the mcsplit case.
            MCSPLIT_FORMAT_FLAG=$(echo "$FORMAT_FLAG" | sed 's/-i/-r/')
            REF_OUT=$(timeout 1000 ./reference/mcsplit/mcsp -b $MCSPLIT_FORMAT_FLAG -q "$A" "$B" 2>/dev/null)
            REF_SIZE=$(echo "$REF_OUT" | grep "Solution size" | awk '{print $3}')
            REF_NODES=$(echo "$REF_OUT" | grep "^Nodes:" | awk '{print $NF}')
            REF_TIME_MS=$(echo "$REF_OUT" | grep "^CPU time" | awk '{print $NF}')
            REF_TIME=$(echo "$REF_TIME_MS" | awk '{printf "%.6f", $1/1000}')
            REF_BIDOMAIN_RAW=$(echo "$REF_OUT" | grep "^BIDOMAINS:" | sed 's/^BIDOMAINS://')
            REF_BIDOMAIN_MEAN=$(echo "$REF_BIDOMAIN_RAW" | cut -d',' -f1)
            REF_BIDOMAIN_MAX=$(echo "$REF_BIDOMAIN_RAW" | cut -d',' -f2)
            ;;
        rl)
            REF_OUT=$(timeout 1000 ./reference/mcsplit-rl/mcsp+RL min_max $FORMAT_FLAG -q "$A" "$B" 2>/dev/null)
            REF_SIZE=$(echo "$REF_OUT" | grep "Solution size" | awk '{print $3}')
            REF_NODES=$(echo "$REF_OUT" | grep "^Nodes:" | awk '{print $NF}')
            REF_TIME_MS=$(echo "$REF_OUT" | grep "^CPU time" | awk '{print $NF}')
            REF_TIME=$(echo "$REF_TIME_MS" | awk '{printf "%.6f", $1/1000}')
            REF_BIDOMAIN_RAW=$(echo "$REF_OUT" | grep "^BIDOMAINS:" | sed 's/^BIDOMAINS://')
            REF_BIDOMAIN_MEAN=$(echo "$REF_BIDOMAIN_RAW" | cut -d',' -f1)
            REF_BIDOMAIN_MAX=$(echo "$REF_BIDOMAIN_RAW" | cut -d',' -f2)
            ;;
        dal)
            REF_OUT=$(timeout 1000 ./reference/mcsplit-dal/mcspDAL min_max $FORMAT_FLAG -q "$A" "$B" 2>/dev/null)
            REF_SIZE=$(echo "$REF_OUT" | grep "^#2:" | awk '{print $2}')
            REF_NODES=$(echo "$REF_OUT" | grep "^#2:" | awk '{print $3}')
            REF_TIME_MS=$(echo "$REF_OUT" | grep "^#2:" | awk '{print $4}')
            REF_TIME=$(echo "$REF_TIME_MS" | awk '{printf "%.6f", $1/1000}')
            REF_BIDOMAIN_RAW=$(echo "$REF_OUT" | grep "^BIDOMAINS:" | sed 's/^BIDOMAINS://')
            REF_BIDOMAIN_MEAN=$(echo "$REF_BIDOMAIN_RAW" | cut -d',' -f1)
            REF_BIDOMAIN_MAX=$(echo "$REF_BIDOMAIN_RAW" | cut -d',' -f2)
            ;;
        ll)
            REF_OUT=$(timeout 1000 ./reference/mcsplit-ll/mcsp+ll min_max $FORMAT_FLAG -q "$A" "$B" 2>/dev/null)
            REF_SIZE=$(echo "$REF_OUT" | grep "Solution size" | awk '{print $3}')
            REF_NODES=$(echo "$REF_OUT" | grep "^Nodes:" | awk '{print $NF}')
            REF_TIME_MS=$(echo "$REF_OUT" | grep "^CPU time" | awk '{print $NF}')
            REF_TIME=$(echo "$REF_TIME_MS" | awk '{printf "%.6f", $1/1000}')
            REF_BIDOMAIN_RAW=$(echo "$REF_OUT" | grep "^BIDOMAINS:" | sed 's/^BIDOMAINS://')
            REF_BIDOMAIN_MEAN=$(echo "$REF_BIDOMAIN_RAW" | cut -d',' -f1)
            REF_BIDOMAIN_MAX=$(echo "$REF_BIDOMAIN_RAW" | cut -d',' -f2)
            ;;
        dsb)
            REF_OUT=$(timeout 1000 ./reference/mcsplit-dsb/bin/run.o min_max $FORMAT_FLAG -q "$A" "$B" 2>/dev/null)
            REF_SIZE=$(echo "$REF_OUT" | awk -F',' '{print $1}' | tr -d ' ')
            REF_TIME=$(echo "$REF_OUT" | awk -F',' '{print $4}' | tr -d ' ')
            REF_NODES=$(echo "$REF_OUT" | awk -F',' '{print $5}' | tr -d ' ')
            REF_BIDOMAIN_RAW=$(echo "$REF_OUT" | grep "^BIDOMAINS:" | sed 's/^BIDOMAINS://')
            REF_BIDOMAIN_MEAN=$(echo "$REF_BIDOMAIN_RAW" | cut -d',' -f1)
            REF_BIDOMAIN_MAX=$(echo "$REF_BIDOMAIN_RAW" | cut -d',' -f2)
            ;;
        rrsplit)
            REF_OUT=$(timeout 1000 ./reference/rrsplit/mcsp min_max $FORMAT_FLAG -q "$A" "$B" 2>/dev/null)
            REF_SIZE=$(echo "$REF_OUT" | grep "^#2:" | awk '{print $2}')
            REF_NODES=$(echo "$REF_OUT" | grep "^#2:" | awk '{print $3}')
            REF_TIME_MS=$(echo "$REF_OUT" | grep "^#2:" | awk '{print $4}')
            REF_TIME=$(echo "$REF_TIME_MS" | awk '{printf "%.6f", $1/1000}')
            REF_BIDOMAIN_RAW=$(echo "$REF_OUT" | grep "^BIDOMAINS:" | sed 's/^BIDOMAINS://')
            REF_BIDOMAIN_MEAN=$(echo "$REF_BIDOMAIN_RAW" | cut -d',' -f1)
            REF_BIDOMAIN_MAX=$(echo "$REF_BIDOMAIN_RAW" | cut -d',' -f2)
            ;;
        symsplit)
            REF_OUT=$(timeout 1000 ./reference/symsplit/bin/run.o min_max $FORMAT_FLAG -q "$A" "$B" 2>/dev/null)
            REF_SIZE=$(echo "$REF_OUT" | awk -F',' '{print $1}' | tr -d ' ')
            REF_TIME=$(echo "$REF_OUT" | awk -F',' '{print $4}' | tr -d ' ')
            REF_NODES=$(echo "$REF_OUT" | awk -F',' '{print $5}' | tr -d ' ')
            REF_BIDOMAIN_RAW=$(echo "$REF_OUT" | grep "^BIDOMAINS:" | sed 's/^BIDOMAINS://')
            REF_BIDOMAIN_MEAN=$(echo "$REF_BIDOMAIN_RAW" | cut -d',' -f1)
            REF_BIDOMAIN_MAX=$(echo "$REF_BIDOMAIN_RAW" | cut -d',' -f2)
            ;;
        *)
            echo "WARNING: No reference binary configured for $ALGO" >&2
            ;;
    esac

	MATCH="UNKNOWN"
	if [ "$UNIFIED_SIZE" == "$REF_SIZE" ]; then
		MATCH="PASS"
	elif [ -z "$UNIFIED_SIZE" ] || [ -z "$REF_SIZE" ]; then
		MATCH="MISSING_DATA"
	else
		MATCH="FAIL"
	fi
fi

# Write output
if [ "$WITH_REF" = "1" ]; then
    echo "$(basename "$A"),$(basename "$B"),$ALGO,$UNIFIED_SIZE,$UNIFIED_EDGES,$UNIFIED_NODES,$UNIFIED_TIME,$UNIFIED_ABORTED,$UNIFIED_NODES_TO_BEST,$UNIFIED_TIME_TO_BEST,$UNIFIED_CUT_BRANCHES,$UNIFIED_BOUND_PRUNED,$UNIFIED_SYM_PRUNED,$UNIFIED_BIDOMAIN_MEAN,$UNIFIED_BIDOMAIN_MAX,$REF_SIZE,$REF_NODES,$REF_TIME,$REF_BIDOMAIN_MEAN,$REF_BIDOMAIN_MAX,$MATCH" > "$OUTFILE"
else
    echo "$(basename "$A"),$(basename "$B"),$ALGO,$UNIFIED_SIZE,$UNIFIED_EDGES,$UNIFIED_NODES,$UNIFIED_TIME,$UNIFIED_ABORTED,$UNIFIED_NODES_TO_BEST,$UNIFIED_TIME_TO_BEST,$UNIFIED_CUT_BRANCHES,$UNIFIED_BOUND_PRUNED,$UNIFIED_SYM_PRUNED,$UNIFIED_BIDOMAIN_MEAN,$UNIFIED_BIDOMAIN_MAX" > "$OUTFILE"
fi
