#!/bin/bash
# submit_lv.sh - submit ONE algorithm for LV
# Usage: bash hpc/cross-dataset-evaluation/scripts/submit_lv.sh mcsplit

ALGO=$1
if [ -z "$ALGO" ]; then
    echo "Usage: bash hpc/cross-dataset-evaluation/scripts/submit_lv.sh <algo>"
    echo "Algos: mcsplit rl ll dal dsb rrsplit symsplit"
    exit 1
fi

mkdir -p hpc/logs hpc/cross-dataset-evaluation/scripts hpc/cross-dataset-evaluation/results

# Generate LV pairs (exclude macOS metadata files)
LIST="hpc/instances_${ALGO}_lv.txt"
FILES=($(ls instances/SIP/LV/g* | grep -v '/\._' | sort -t'g' -k2 -n))
> "$LIST"
for ((i=0; i<${#FILES[@]}; i++)); do
    for ((j=i+1; j<${#FILES[@]}; j++)); do
        echo "${FILES[$i]} ${FILES[$j]}"
    done
done >> "$LIST"
echo "LV graphs: ${#FILES[@]}"
echo "LV pairs: $(wc -l < $LIST)"

NJOBS=$(($(wc -l < "$LIST") - 1))
SCRIPT="hpc/cross-dataset-evaluation/scripts/job_${ALGO}_lv.slurm"
sed "s/{ALGO}/$ALGO/g; s/{NJOBS}/$NJOBS/g; s/{DATASET}/lv/g" \
    hpc/cross-dataset-evaluation/scripts/job_array.slurm > "$SCRIPT"
sed -i 's/\r$//' "$SCRIPT"

if grep -q "{DATASET}" "$SCRIPT"; then
    echo "ERROR: {DATASET} was not replaced in $SCRIPT"
    exit 1
fi

echo "Submitting $ALGO LV: $((NJOBS+1)) jobs"
sbatch "$SCRIPT"