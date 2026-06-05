# Reference binaries
REF_MCSPLIT="./reference/mcsplit/mcsp"
REF_DAL="./reference/mcsplit-dal/mcspDAL"
REF_DSB="./reference/mcsplit-dsb/bin/run.o"
REF_LL="./reference/mcsplit-ll/mcsp+ll"
REF_RL="./reference/mcsplit-rl/mcsp+RL"
REF_RRSPLIT="./reference/rrsplit/mcsp"
REF_SYMSPLIT="./reference/symsplit/bin/run.o"

# Unified binary
UNIFIED="./solvers/bin/mcsp"

# Output goes to results/
OUTPUT="results/validation_results_${ALGO}.csv"

# Parse -f flag
SINGLE_FILE=""
if [ "$1" == "-f" ]; then
    ALGO="$2"
    SINGLE_FILE="$3"
else
    ALGO="${1:-mcsplit}"
    INSTANCE_DIR="${2:-instances/MIVIA/mcs10/bvg/b03}"
    NUM="${3:-1}"
    PATTERN="${4:-*.A00}"
fi

# Select reference binary and command based on algorithm
case "$ALGO" in
    mcsplit)
        REF_BIN="$REF_MCSPLIT"
        REF_CMD() { $REF_BIN -q "$1" "$2" 2>/dev/null; }
        parse_ref() {
            ref_size=$(echo "$1" | grep "Solution size" | awk '{print $3}')
            ref_edges=$(echo "$1" | grep "Number of edges" | awk '{print $6}')
            ref_nodes=$(echo "$1" | grep "Nodes:" | awk '{print $2}' | tr -d ',')
        }
        ;;
    dal)
        REF_BIN="$REF_DAL"
        REF_CMD() { $REF_BIN min_max -q "$1" "$2" 2>/dev/null; }
        parse_ref() {
            ref_size=$(echo "$1" | grep "^#2:" | awk '{print $2}')
            ref_edges="N/A"
            ref_nodes=$(echo "$1" | grep "^#2:" | awk '{print $3}')
        }
        ;;
    dsb)
        REF_BIN="$REF_DSB"
        REF_CMD() { $REF_BIN min_max -q "$1" "$2" 2>/dev/null; }
        parse_ref() {
            ref_size=$(echo "$1" | awk -F',' '{print $1}')
            ref_edges="N/A"
            ref_nodes=$(echo "$1" | awk -F',' '{print $5}')
        }
        ;;
    ll)
        REF_BIN="$REF_LL"
        REF_CMD() { $REF_BIN min_max -q "$1" "$2" 2>/dev/null; }
        parse_ref() {
            ref_size=$(echo "$1" | grep "Solution size" | awk '{print $3}')
            ref_edges="N/A"
            ref_nodes=$(echo "$1" | grep "Nodes:" | awk '{print $2}' | tr -d ',')
        }
        ;;
    rl)
        REF_BIN="$REF_RL"
        REF_CMD() { $REF_BIN min_max -q "$1" "$2" 2>/dev/null; }
        parse_ref() {
            ref_size=$(echo "$1" | grep "Solution size" | awk '{print $3}')
            ref_edges="N/A"
            ref_nodes=$(echo "$1" | grep "Nodes:" | awk '{print $2}' | tr -d ',')
        }
        ;;
    rrsplit)
        REF_BIN="$REF_RRSPLIT"
        REF_CMD() { $REF_BIN min_max -q "$1" "$2" 2>/dev/null; }
        parse_ref() {
            ref_size=$(echo "$1" | grep "^#2:" | awk '{print $2}')
            ref_edges="N/A"
            ref_nodes=$(echo "$1" | grep "^#2:" | awk '{print $3}')
        }
        ;;
    symsplit)
        REF_BIN="$REF_SYMSPLIT"
        REF_CMD() { $REF_BIN min_max -q "$1" "$2" 2>/dev/null; }
        parse_ref() {
            ref_size=$(echo "$1" | awk -F',' '{print $1}')
            ref_edges="N/A"
            ref_nodes=$(echo "$1" | awk -F',' '{print $5}')
        }
        ;;
    unified)
        REF_BIN="$UNIFIED"
        REF_CMD() { $REF_BIN -q "$1" "$2" 2>/dev/null; }
        parse_ref() {
            ref_size=$(echo "$1" | awk '{print $1}')
            ref_edges=$(echo "$1" | awk '{print $2}')
            ref_nodes=$(echo "$1" | awk '{print $3}')
        }
        ;;
    *)
        echo "Unknown algorithm: $ALGO"
        echo "Choose from: mcsplit, dal, dsb, ll, rl, rrsplit, symsplit, unified"
        exit 1
        ;;
esac

echo "instance,ref_size,ref_edges,ref_nodes,unified_size,unified_edges,unified_nodes,unified_time,unified_aborted,unified_root_ub,unified_nodes_to_best,unified_time_to_best,unified_cut_branches,unified_bound_pruned,unified_sym_pruned,unified_conflicts,match" > $OUTPUT

count=0
pass=0
fail=0

echo "=========================================="
echo "Validating unified vs reference: $ALGO"
echo "Reference binary:  $REF_BIN"
echo "Unified binary:    $UNIFIED"
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
    echo "  -> running reference ($ALGO)..."
    ref_out=$(REF_CMD "$A" "$B")
    parse_ref "$ref_out"
    echo "  -> reference done: size=$ref_size edges=$ref_edges nodes=$ref_nodes"

    echo "  -> running unified..."
    unified_out=$($UNIFIED -q "$A" "$B" 2>/dev/null)
    unified_size=$(echo "$unified_out"      | awk '{print $1}')
    unified_edges=$(echo "$unified_out"     | awk '{print $2}')
    unified_nodes=$(echo "$unified_out"     | awk '{print $3}')
    unified_time=$(echo "$unified_out"      | awk '{print $4}')
    unified_aborted=$(echo "$unified_out"   | awk '{print $5}')
    unified_root_ub=$(echo "$unified_out"   | awk '{print $6}')
    unified_n2best=$(echo "$unified_out"    | awk '{print $7}')
    unified_t2best=$(echo "$unified_out"    | awk '{print $8}')
    unified_cut=$(echo "$unified_out"       | awk '{print $9}')
    unified_bound=$(echo "$unified_out"     | awk '{print $10}')
    unified_sym=$(echo "$unified_out"       | awk '{print $11}')
    unified_conflicts=$(echo "$unified_out" | awk '{print $12}')
    echo "  -> unified done: size=$unified_size edges=$unified_edges nodes=$unified_nodes time=${unified_time}s"

    # Compare — skip edge comparison if reference doesn't report edges
    if [ "$ref_size" == "$unified_size" ] && [ "$ref_nodes" == "$unified_nodes" ]; then
        if [ "$ref_edges" == "N/A" ] || [ "$ref_edges" == "$unified_edges" ]; then
            match="PASS"
            pass=$((pass + 1))
            echo "  -> PASS (size and nodes match)"
        else
            match="PASS_SIZE_NODES"
            pass=$((pass + 1))
            echo "  -> PASS (size and nodes match) WARNING: edges differ ref=$ref_edges unified=$unified_edges"
        fi
    elif [ "$ref_size" == "$unified_size" ]; then
        match="PASS_SIZE_ONLY"
        pass=$((pass + 1))
        echo "  -> PASS (size matches) WARNING: nodes differ ref=$ref_nodes unified=$unified_nodes"
    else
        match="FAIL"
        fail=$((fail + 1))
        echo "  -> FAIL *** MISMATCH ref=$ref_size unified=$unified_size ***"
    fi

    echo "$instance,$ref_size,$ref_edges,$ref_nodes,$unified_size,$unified_edges,$unified_nodes,$unified_time,$unified_aborted,$unified_root_ub,$unified_n2best,$unified_t2best,$unified_cut,$unified_bound,$unified_sym,$unified_conflicts,$match" >> $OUTPUT
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