#!/bin/bash
# validate_all_agree.sh
ALGOS="mcsplit dal ll dsb rrsplit symsplit"
INSTANCE_DIR=${1:-instances/MIVIA/mcs10/rand/r01}
SIZES="${2:-s20 s25 s30 s35 s40}"
TIMEOUT=${3:-5}
MAX_INSTANCES=${4:-0}   # 0 = no limit, applies per size
GTR_FILE="${INSTANCE_DIR}.gtr"

mkdir -p results
OUTFILE="results/cross_validation.csv"
echo "instance,ground_truth,$(echo $ALGOS | tr ' ' ',')" > $OUTFILE

declare -A GT
HAS_GTR=false
if [ -f "$GTR_FILE" ]; then
    HAS_GTR=true
    while read -r name size; do
        GT["$name"]="$size"
    done < "$GTR_FILE"
    echo "Loaded ground truth from $GTR_FILE"
else
    echo "No ground truth file found at $GTR_FILE — running cross-solver agreement only"
fi

# Count total files across all requested sizes
TOTAL_FILES=0
for SZ in $SIZES; do
    n=$(ls $INSTANCE_DIR/*${SZ}*.A* 2>/dev/null | wc -l)
    if [ "$MAX_INSTANCES" -gt 0 ] && [ "$n" -gt "$MAX_INSTANCES" ]; then
        n=$MAX_INSTANCES
    fi
    TOTAL_FILES=$((TOTAL_FILES + n))
done
echo "Will check $TOTAL_FILES instances across sizes: $SIZES"
echo ""

mismatch_count=0
gt_mismatch_count=0
total_count=0

for SZ in $SIZES; do
    size_count=0
    for A in $INSTANCE_DIR/*${SZ}*.A*; do
        [ -f "$A" ] || continue
        B="${A/.A/.B}"
        [ -f "$B" ] || continue

        if [ "$MAX_INSTANCES" -gt 0 ] && [ "$size_count" -ge "$MAX_INSTANCES" ]; then
            break
        fi

        instance=$(basename "$A")
        total_count=$((total_count + 1))
        size_count=$((size_count + 1))

        printf "[%d/%d] %s ... " "$total_count" "$TOTAL_FILES" "$instance"

        gt_size="N/A"
        if [ "$HAS_GTR" = true ]; then
            gt_size="${GT[$instance]:-N/A}"
        fi

        sizes=""
        for algo in $ALGOS; do
            size=$(./solvers/bin/mcsp -A $algo -t $TIMEOUT -q "$A" "$B" 2>/dev/null | awk '{print $1}')
            sizes="$sizes,$size"
        done

        echo "$instance,$gt_size$sizes" >> $OUTFILE

        unique=$(echo "$sizes" | tr ',' '\n' | grep -v '^$' | sort -u | wc -l)
        if [ "$unique" -gt 1 ]; then
            echo "MISMATCH: $sizes"
            mismatch_count=$((mismatch_count + 1))
        else
            echo "OK (size=$(echo $sizes | cut -d',' -f2))"
        fi

        if [ "$gt_size" != "N/A" ]; then
            first_size=$(echo "$sizes" | cut -d',' -f2)
            if [ "$first_size" != "$gt_size" ]; then
                echo "  GT MISMATCH: ground_truth=$gt_size got=$first_size"
                gt_mismatch_count=$((gt_mismatch_count + 1))
            fi
        fi
    done
done

echo ""
echo "=========================================="
echo "Checked $total_count instances"
echo "Cross-solver mismatches: $mismatch_count"
echo "Ground-truth mismatches: $gt_mismatch_count"
echo "Results saved to $OUTFILE"
echo "=========================================="