#!/bin/bash
# Usage: bash hpc/test_references.sh <file_A> <file_B>

A=$1
B=$2

if [ -z "$A" ] || [ -z "$B" ]; then
    echo "Usage: bash hpc/test_references.sh <file_A> <file_B>"
    exit 1
fi

echo "Instance: $A $B"
echo ""

echo "============================================================"
echo "mcsplit"
echo "============================================================"
./reference/mcsplit/mcsp -q "$A" "$B" 2>/dev/null

echo ""
echo "============================================================"
echo "rl"
echo "============================================================"
./reference/mcsplit-rl/mcsp+RL min_max -q "$A" "$B" 2>/dev/null

echo ""
echo "============================================================"
echo "dal"
echo "============================================================"
./reference/mcsplit-dal/mcspDAL min_max -q "$A" "$B" 2>/dev/null

echo ""
echo "============================================================"
echo "ll"
echo "============================================================"
./reference/mcsplit-ll/mcsp+ll min_max -q "$A" "$B" 2>/dev/null

echo ""
echo "============================================================"
echo "dsb"
echo "============================================================"
./reference/mcsplit-dsb/bin/run.o min_max -q "$A" "$B" 2>/dev/null

echo ""
echo "============================================================"
echo "rrsplit"
echo "============================================================"
./reference/rrsplit/mcsp min_max -q "$A" "$B" 2>/dev/null

echo ""
echo "============================================================"
echo "symsplit"
echo "============================================================"
./reference/symsplit/bin/run.o min_max -q "$A" "$B" 2>/dev/null