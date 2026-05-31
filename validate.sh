#!/bin/bash
# validate.sh - Compare unified McSplit against reference implementation
# Usage: ./validate.sh [instance_dir] [num_instances] [pattern]
#        ./validate.sh -f instances/MIVIA/mcs10/bvg/b03/mcs10_b03_s40.A00

UNIFIED="./mcsp"
REFERENCE="./reference/mcsplit/code/james-c/mcsp"
INSTANCE_DIR="${1:-instances/MIVIA/mcs10/bvg/b03}"
NUM="${2:-1}"
PATTERN="${3:-*.A00}"
OUTPUT="validation_results.csv"
SINGLE_FILE=""

# Check for -f flag
if [ "$1" == "-f" ]; then
    SINGLE_FILE="$2"
fi

echo "instance,ref_size,unified_size,ref_nodes,unified_nodes,unified_time,match" > $OUTPUT

count=0
pass=0
fail=0

echo "=========================================="
echo "Validating McSplit unified vs reference"
echo "Unified binary:    $UNIFIED"
echo "Reference binary:  $REFERENCE"
echo "=========================================="
echo ""

run_instance() {
    local A=$1
    local B="${A/.A00/.B00}"
    if [ ! -f "$B" ]; then
        echo "  No matching .B00 for $A"
        return
    fi

    local instance=$(basename $A .A00)
    count=$((count + 1))

    echo "[$count/$NUM] $instance"
    echo "  -> running reference..."
    ref_out=$($REFERENCE -q $A $B 2>/dev/null)
    ref_size=$(echo "$ref_out" | grep "Solution size" | awk '{print $3}')
    ref_nodes=$(echo "$ref_out" | grep "Nodes:" | awk '{print $2}' | tr -d ',')
    echo "  -> reference done: size=$ref_size nodes=$ref_nodes"

    echo "  -> running unified..."
    unified_out=$($UNIFIED -q $A $B 2>/dev/null)
    unified_size=$(echo "$unified_out" | awk '{print $1}')
    unified_nodes=$(echo "$unified_out" | awk '{print $2}')
    unified_time=$(echo "$unified_out" | awk '{print $3}')
    echo "  -> unified done: size=$unified_size nodes=$unified_nodes time=${unified_time}s"

    if [ "$ref_size" == "$unified_size" ]; then
        match="PASS"
        pass=$((pass + 1))
        echo "  -> PASS"
    else
        match="FAIL"
        fail=$((fail + 1))
        echo "  -> FAIL *** MISMATCH ref=$ref_size unified=$unified_size ***"
    fi

    echo "$instance,$ref_size,$unified_size,$ref_nodes,$unified_nodes,$unified_time,$match" >> $OUTPUT
    echo ""
}

if [ -n "$SINGLE_FILE" ]; then
    echo "Single file: $SINGLE_FILE"
    echo ""
    run_instance "$SINGLE_FILE"
else
    echo "Instance dir: $INSTANCE_DIR"
    echo "Pattern:      $PATTERN"
    echo "Running up to $NUM instances"
    echo ""
    for A in $INSTANCE_DIR/$PATTERN; do
        if [ ! -f "$A" ]; then
            echo "No files matching pattern $PATTERN in $INSTANCE_DIR"
            break
        fi
        run_instance "$A"
        if [ $count -ge $NUM ]; then break; fi
    done
fi

echo "=========================================="
echo "Results: $pass/$count passed, $fail failed"
echo "=========================================="
echo "Saved to $OUTPUT"