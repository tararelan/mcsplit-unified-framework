#!/bin/bash
ALGO=$1
A=$2
B=$3

OUTDIR="results/hpc/${ALGO}"
mkdir -p "$OUTDIR"

OUTFILE="${OUTDIR}/$(basename "$A").$(basename "$B").csv"

if [ -f "$OUTFILE" ]; then
    exit 0  # already done, skip
fi

RESULT=$(./solvers/bin/mcsp -A "$ALGO" -t 1000 -q "$A" "$B" 2>/dev/null)
echo "$(basename "$A"),$(basename "$B"),$RESULT" > "$OUTFILE"