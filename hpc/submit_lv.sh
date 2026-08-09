#!/bin/bash
# submit_lv.sh - submit one algorithm for the LV dataset
# Usage: bash hpc/submit_lv.sh <algo> [--with-ref]
# --with-ref: also run reference binary and compare

ALGO=$1
WITH_REF=0
if [ "$2" = "--with-ref" ]; then
    WITH_REF=1
fi

if [ -z "$ALGO" ]; then
    echo "Usage: bash hpc/submit_lv.sh <algo> [--with-ref]"
    echo "Algorithms: mcsplit rl ll ll_lsm ll_lum dal dal_dal dal_rl dsb dsb_always dsb_never rrsplit rrsplit_noveq rrsplit_nomax rrsplit_nobound symsplit symsplit_valonly symsplit_varonly"
    exit 1
fi

mkdir -p hpc/logs hpc/results

LIST="hpc/instances_${ALGO}_lv.txt"
FILES=($(ls instances/SIP/LV/g* | grep -v '/\._' | sort -t'g' -k2 -n))
> "$LIST"
for ((i=0; i<${#FILES[@]}; i++)); do
    for ((j=i+1; j<${#FILES[@]}; j++)); do
        echo "${FILES[$i]} ${FILES[$j]}"
    done
done >> "$LIST"
echo "LV pairs: $(wc -l < $LIST)"

NJOBS=$(($(wc -l < "$LIST") - 1))

if [ "$WITH_REF" = "1" ]; then
    SCRIPT="hpc/job_${ALGO}_lv_ref.slurm"
else
    SCRIPT="hpc/job_${ALGO}_lv.slurm"
fi

sed "s/{ALGO}/$ALGO/g; s/{NJOBS}/$NJOBS/g; s/{DATASET}/lv/g; s/{WITH_REF}/$WITH_REF/g; s#{LIST}#$LIST#g" \
    hpc/job_array.slurm > "$SCRIPT"
sed -i 's/\r$//' "$SCRIPT"

if grep -qE "\{ALGO\}|\{NJOBS\}|\{DATASET\}|\{LIST\}|\{WITH_REF\}" "$SCRIPT"; then
    echo "ERROR: placeholder not replaced in $SCRIPT"
    exit 1
fi

if [ "$WITH_REF" = "1" ]; then
    echo "Submitting $ALGO LV (with reference comparison): $((NJOBS+1)) jobs"
else
    echo "Submitting $ALGO LV: $((NJOBS+1)) jobs"
fi
sbatch "$SCRIPT"