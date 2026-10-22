#!/bin/bash
# submit_bi.sh - submit one algorithm for the BI dataset
# Usage: bash hpc/submit_bi.sh <algo> [--with-ref] [--swap]
# --with-ref: also run reference binary and compare
# --swap:     uniform preprocessing - smaller graph first for EVERY algorithm
#             (solver -s flag); results go to hpc/results/..._swap/

ALGO=$1
WITH_REF=0
SWAP=0
for opt in "${@:2}"; do
    case "$opt" in
        --with-ref) WITH_REF=1 ;;
        --swap)     SWAP=1 ;;
        *) echo "Unknown option: $opt"; exit 1 ;;
    esac
done

if [ -z "$ALGO" ]; then
    echo "Usage: bash hpc/submit_bi.sh <algo> [--with-ref]"
    echo "Algorithms: mcsplit rl ll ll_lsm ll_lum dal dal_dal dal_rl dsb dsb_always dsb_never rrsplit rrsplit_noveq rrsplit_nomax rrsplit_nobound symsplit symsplit_valonly symsplit_varonly"
    exit 1
fi

mkdir -p hpc/logs hpc/results

LIST="hpc/instances_${ALGO}_bi.txt"
FILES=($(ls instances/SIP/biochemicalReactions/*.txt | sort))
> "$LIST"
for ((i=0; i<${#FILES[@]}; i++)); do
    for ((j=i+1; j<${#FILES[@]}; j++)); do
        echo "${FILES[$i]} ${FILES[$j]}"
    done
done >> "$LIST"
echo "BI pairs: $(wc -l < $LIST)"

NJOBS=$(($(wc -l < "$LIST") - 1))

if [ "$WITH_REF" = "1" ]; then
    SCRIPT="hpc/job_${ALGO}_bi_ref.slurm"
else
    SCRIPT="hpc/job_${ALGO}_bi.slurm"
fi

# separate generated script per mode so swap and faithful jobs don't overwrite each other
if [ "$SWAP" = "1" ]; then
    SCRIPT="${SCRIPT%.slurm}_swap.slurm"
fi

sed "s/{ALGO}/$ALGO/g; s/{NJOBS}/$NJOBS/g; s/{DATASET}/bi/g; s/{WITH_REF}/$WITH_REF/g; s#{LIST}#$LIST#g" \
    hpc/job_array.slurm | sed "s/{SWAP}/$SWAP/g" > "$SCRIPT"
sed -i 's/\r$//' "$SCRIPT"

if grep -qE "\{ALGO\}|\{NJOBS\}|\{DATASET\}|\{LIST\}|\{WITH_REF\}|\{SWAP\}" "$SCRIPT"; then
    echo "ERROR: placeholder not replaced in $SCRIPT"
    exit 1
fi

if [ "$WITH_REF" = "1" ]; then
    echo "Submitting $ALGO BI (with reference comparison): $((NJOBS+1)) jobs"
else
    echo "Submitting $ALGO BI: $((NJOBS+1)) jobs"
fi
[ "$SWAP" = "1" ] && echo "  (uniform preprocessing: --swap)"
sbatch "$SCRIPT"