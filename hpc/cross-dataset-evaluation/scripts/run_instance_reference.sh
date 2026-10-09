#!/bin/bash
ALGO=$1
A=$2
B=$3
DATASET=$4  # bi or lv

OUTDIR="hpc/cross-dataset-evaluation/results/${ALGO}_${DATASET}_reference"
mkdir -p "$OUTDIR"

OUTFILE="${OUTDIR}/$(basename "$A").$(basename "$B").csv"

if [ -f "$OUTFILE" ]; then
    exit 0
fi

if [ "$DATASET" = "bi" ]; then
    FORMAT_FLAG="-l -i"
elif [ "$DATASET" = "lv" ]; then
    FORMAT_FLAG="-l"
fi

# Reference binaries: mcsplit's directed flag is -r, all others use -i
REF_DIRECTED_FLAG=""
if [ "$DATASET" = "bi" ]; then
    if [ "$ALGO" = "mcsplit" ]; then
        REF_DIRECTED_FLAG="-r"
    else
        REF_DIRECTED_FLAG="-i"
    fi
fi

# UNIFIED_RESULT=$(./solvers/bin/mcsp -A "$ALGO" -t 1000 $FORMAT_FLAG -q "$A" "$B" 2>/dev/null)

REF_SIZE=""
REF_NODES=""
REF_TIME=""
REF_BIDOMAIN_MEAN=""
REF_BIDOMAIN_MAX=""

case "$ALGO" in
    mcsplit)
        REF_OUT=$(timeout 1000 ./reference/mcsplit/mcsp $REF_DIRECTED_FLAG -l -q "$A" "$B" 2>/dev/null)
        REF_SIZE=$(echo "$REF_OUT" | grep "Solution size" | awk '{print $3}')
        REF_NODES=$(echo "$REF_OUT" | grep "^Nodes:" | awk '{print $NF}')
        REF_TIME_MS=$(echo "$REF_OUT" | grep "^CPU time" | awk '{print $NF}')
        REF_TIME=$(echo "$REF_TIME_MS" | awk '{printf "%.6f", $1/1000}')
        REF_BIDOMAIN_RAW=$(echo "$REF_OUT" | grep "^BIDOMAINS:" | sed 's/^BIDOMAINS://')
        REF_BIDOMAIN_MEAN=$(echo "$REF_BIDOMAIN_RAW" | cut -d',' -f1)
        REF_BIDOMAIN_MAX=$(echo "$REF_BIDOMAIN_RAW" | cut -d',' -f2)
        ;;
    rl)
        REF_OUT=$(timeout 1000 ./reference/mcsplit-rl/mcsp+RL min_max $REF_DIRECTED_FLAG -l -q "$A" "$B" 2>/dev/null)
        REF_SIZE=$(echo "$REF_OUT" | grep "Solution size" | awk '{print $3}')
        REF_NODES=$(echo "$REF_OUT" | grep "^Nodes:" | awk '{print $NF}')
        REF_TIME_MS=$(echo "$REF_OUT" | grep "^CPU time" | awk '{print $NF}')
        REF_TIME=$(echo "$REF_TIME_MS" | awk '{printf "%.6f", $1/1000}')
        REF_BIDOMAIN_RAW=$(echo "$REF_OUT" | grep "^BIDOMAINS:" | sed 's/^BIDOMAINS://')
        REF_BIDOMAIN_MEAN=$(echo "$REF_BIDOMAIN_RAW" | cut -d',' -f1)
        REF_BIDOMAIN_MAX=$(echo "$REF_BIDOMAIN_RAW" | cut -d',' -f2)
        ;;
    dal)
        REF_OUT=$(timeout 1000 ./reference/mcsplit-dal/mcspDAL min_max $REF_DIRECTED_FLAG -l -l -q "$A" "$B" 2>/dev/null)
        REF_SIZE=$(echo "$REF_OUT" | grep "^#2:" | awk '{print $2}')
        REF_NODES=$(echo "$REF_OUT" | grep "^#2:" | awk '{print $3}')
        REF_TIME_MS=$(echo "$REF_OUT" | grep "^#2:" | awk '{print $4}')
        REF_TIME=$(echo "$REF_TIME_MS" | awk '{printf "%.6f", $1/1000}')
        REF_BIDOMAIN_RAW=$(echo "$REF_OUT" | grep "^BIDOMAINS:" | sed 's/^BIDOMAINS://')
        REF_BIDOMAIN_MEAN=$(echo "$REF_BIDOMAIN_RAW" | cut -d',' -f1)
        REF_BIDOMAIN_MAX=$(echo "$REF_BIDOMAIN_RAW" | cut -d',' -f2)
        ;;
    ll)
        REF_OUT=$(timeout 1000 ./reference/mcsplit-ll/mcsp+ll min_max $REF_DIRECTED_FLAG -l -q "$A" "$B" 2>/dev/null)
        REF_SIZE=$(echo "$REF_OUT" | grep "Solution size" | awk '{print $3}')
        REF_NODES=$(echo "$REF_OUT" | grep "^Nodes:" | awk '{print $NF}')
        REF_TIME_MS=$(echo "$REF_OUT" | grep "^CPU time" | awk '{print $NF}')
        REF_TIME=$(echo "$REF_TIME_MS" | awk '{printf "%.6f", $1/1000}')
        REF_BIDOMAIN_RAW=$(echo "$REF_OUT" | grep "^BIDOMAINS:" | sed 's/^BIDOMAINS://')
        REF_BIDOMAIN_MEAN=$(echo "$REF_BIDOMAIN_RAW" | cut -d',' -f1)
        REF_BIDOMAIN_MAX=$(echo "$REF_BIDOMAIN_RAW" | cut -d',' -f2)
        ;;
    dsb)
        REF_OUT=$(timeout 1000 ./reference/mcsplit-dsb/bin/run.o min_max $REF_DIRECTED_FLAG -l -q "$A" "$B" 2>/dev/null)
        REF_SIZE=$(echo "$REF_OUT" | awk -F',' '{print $1}' | tr -d ' ')
        REF_TIME=$(echo "$REF_OUT" | awk -F',' '{print $4}' | tr -d ' ')
        REF_NODES=$(echo "$REF_OUT" | awk -F',' '{print $5}' | tr -d ' ')
        REF_BIDOMAIN_RAW=$(echo "$REF_OUT" | grep "^BIDOMAINS:" | sed 's/^BIDOMAINS://')
        REF_BIDOMAIN_MEAN=$(echo "$REF_BIDOMAIN_RAW" | cut -d',' -f1)
        REF_BIDOMAIN_MAX=$(echo "$REF_BIDOMAIN_RAW" | cut -d',' -f2)
        ;;
    rrsplit)
        REF_OUT=$(timeout 1000 ./reference/rrsplit/mcsp min_max $REF_DIRECTED_FLAG -l -q "$A" "$B" 2>/dev/null)
        REF_SIZE=$(echo "$REF_OUT" | grep "^#2:" | awk '{print $2}')
        REF_NODES=$(echo "$REF_OUT" | grep "^#2:" | awk '{print $3}')
        REF_TIME_MS=$(echo "$REF_OUT" | grep "^#2:" | awk '{print $4}')
        REF_TIME=$(echo "$REF_TIME_MS" | awk '{printf "%.6f", $1/1000}')
        REF_BIDOMAIN_RAW=$(echo "$REF_OUT" | grep "^BIDOMAINS:" | sed 's/^BIDOMAINS://')
        REF_BIDOMAIN_MEAN=$(echo "$REF_BIDOMAIN_RAW" | cut -d',' -f1)
        REF_BIDOMAIN_MAX=$(echo "$REF_BIDOMAIN_RAW" | cut -d',' -f2)
        ;;
    symsplit)
        REF_OUT=$(timeout 1000 ./reference/symsplit/bin/run.o min_max $REF_DIRECTED_FLAG -l -q "$A" "$B" 2>/dev/null)
        REF_SIZE=$(echo "$REF_OUT" | awk -F',' '{print $1}' | tr -d ' ')
        REF_TIME=$(echo "$REF_OUT" | awk -F',' '{print $4}' | tr -d ' ')
        REF_NODES=$(echo "$REF_OUT" | awk -F',' '{print $5}' | tr -d ' ')
        REF_BIDOMAIN_RAW=$(echo "$REF_OUT" | grep "^BIDOMAINS:" | sed 's/^BIDOMAINS://')
        REF_BIDOMAIN_MEAN=$(echo "$REF_BIDOMAIN_RAW" | cut -d',' -f1)
        REF_BIDOMAIN_MAX=$(echo "$REF_BIDOMAIN_RAW" | cut -d',' -f2)
        ;;
esac

echo "$(basename "$A"),$(basename "$B"),$ALGO,$REF_SIZE,$REF_NODES,$REF_TIME,$REF_BIDOMAIN_MEAN,$REF_BIDOMAIN_MAX" > "$OUTFILE"