#!/bin/bash
# submit_bi.sh - submit one algorithm for the BI dataset
# Usage: bash hpc/submit_bi.sh <algo> [--with-ref]
# --with-ref: also run reference binary and compare

ALGO=$1
WITH_REF=0
if [ "$2" = "--with-ref" ]; then
    WITH_REF=1
fi

if [ -z "$ALGO" ]; then
    echo "Usage: bash hpc/submit_bi.sh <algo> [--with-ref]"
    echo "Algorithms: mcsplit rl ll ll_lsm ll_lum dal dal_dal dal_rl dsb dsb_always dsb_never rrsplit rrsplit_noveq rrsplit_nomax rrsplit_nobound symsplit symsplit_valonly symsplit_varonly"
    exit 1
fi

mkdir -p hpc/logs hpc/scripts hpc/results

# Generate BI pairs
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
    SCRIPT="hpc/scripts/job_${ALGO}_bi_ref.slurm"
    sed "s/{ALGO}/$ALGO/g; s/{NJOBS}/$NJOBS/g; s/{DATASET}/bi/g" \
        hpc/scripts/job_array_ref.slurm > "$SCRIPT"
    sed -i 's/\r$//' "$SCRIPT"
    echo "Submitting $ALGO BI (with reference comparison): $((NJOBS+1)) jobs"
else
    SCRIPT="hpc/scripts/job_${ALGO}_bi.slurm"
    sed "s/{ALGO}/$ALGO/g; s/{NJOBS}/$NJOBS/g; s/{DATASET}/bi/g" \
        hpc/scripts/job_array.slurm > "$SCRIPT"
    sed -i 's/\r$//' "$SCRIPT"
    echo "Submitting $ALGO BI: $((NJOBS+1)) jobs"
fi

if grep -q "{DATASET}" "$SCRIPT"; then
    echo "ERROR: {DATASET} was not replaced in $SCRIPT"
    exit 1
fi

sbatch "$SCRIPT"