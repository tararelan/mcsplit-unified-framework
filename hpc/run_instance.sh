#!/bin/bash
ALGO=$1
A=$2
B=$3

OUTDIR="results/hpc/${ALGO}"
mkdir -p $OUTDIR

OUTFILE="${OUTDIR}/$(basename $A .A00)__$(basename $B .B00).csv"

if [ -f "$OUTFILE" ]; then
    exit 0  # already done, skip
fi

RESULT=$(./solvers/bin/mcsp -A $ALGO -t 600 -q "$A" "$B" 2>/dev/null)
echo "$(basename $A),$(basename $B),$RESULT" > $OUTFILE