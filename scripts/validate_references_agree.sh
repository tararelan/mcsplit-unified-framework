#!/bin/bash
# validate_references_agree.sh
# Cross-checks the REFERENCE binaries against each other (not your unified implementation)

REF_MCSPLIT="./reference/mcsplit/mcsp"
REF_DAL="./reference/mcsplit-dal/mcspDAL"
REF_DSB="./reference/mcsplit-dsb/bin/run.o"
REF_LL="./reference/mcsplit-ll/mcsp+ll"
REF_RL="./reference/mcsplit-rl/mcsp+RL"
REF_RRSPLIT="./reference/rrsplit/mcsp"
REF_SYMSPLIT="./reference/symsplit/bin/run.o"

INSTANCE_DIR=${1:-instances/MIVIA/mcs10/rand/r01}
SIZES="${2:-s20 s25 s30 s35 s40}"
TIMEOUT=${3:-5}
MAX_INSTANCES=${4:-0}   # 0 = no limit, applies per size

mkdir -p results
OUTFILE="results/reference_cross_validation.csv"
echo "instance,mcsplit,dal,dsb,ll,rl,rrsplit,symsplit" > $OUTFILE

# Each algorithm has its own command + parser since reference binaries
# all have different CLIs and output formats.
get_size() {
    local algo=$1
    local A=$2
    local B=$3
    local out

    case "$algo" in
        mcsplit)
            out=$($REF_MCSPLIT -q "$A" "$B" 2>/dev/null)
            echo "$out" | grep "Solution size" | awk '{print $3}'
            ;;
        dal)
            out=$(timeout $TIMEOUT $REF_DAL min_max -q "$A" "$B" 2>/dev/null)
            echo "$out" | grep "^#2:" | awk '{print $2}'
            ;;
        dsb)
            out=$(timeout $TIMEOUT $REF_DSB min_max -q "$A" "$B" 2>/dev/null)
            echo "$out" | awk -F',' '{print $1}' | tr -d ' '
            ;;
        ll)
            out=$(timeout $TIMEOUT $REF_LL min_max -q "$A" "$B" 2>/dev/null)
            echo "$out" | grep "Solution size" | awk '{print $3}'
            ;;
        rl)
            out=$(timeout $TIMEOUT $REF_RL min_max -q "$A" "$B" 2>/dev/null)
            echo "$out" | grep "Solution size" | awk '{print $3}'
            ;;
        rrsplit)
            out=$(timeout $TIMEOUT $REF_RRSPLIT min_max -q "$A" "$B" 2>/dev/null)
            echo "$out" | grep "^#2:" | awk '{print $2}'
            ;;
        symsplit)
            out=$(timeout $TIMEOUT $REF_SYMSPLIT min_max -q "$A" "$B" 2>/dev/null)
            echo "$out" | awk -F',' '{print $1}'
            ;;
    esac
}

ALGOS="mcsplit dal dsb ll rl rrsplit symsplit"

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

        sizes=""
        for algo in $ALGOS; do
            size=$(get_size "$algo" "$A" "$B")
            sizes="$sizes,$size"
        done

        echo "$instance$sizes" >> $OUTFILE

        unique=$(echo "$sizes" | tr ',' '\n' | grep -v '^$' | sort -u | wc -l)
        if [ "$unique" -gt 1 ]; then
            echo "MISMATCH: $sizes"
            mismatch_count=$((mismatch_count + 1))
        else
            echo "OK (size=$(echo $sizes | cut -d',' -f2))"
        fi
    done
done

echo ""
echo "=========================================="
echo "Checked $total_count instances"
echo "Reference cross-solver mismatches: $mismatch_count"
echo "Results saved to $OUTFILE"
echo "=========================================="