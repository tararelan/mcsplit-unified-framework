#!/bin/bash
# submit_bi.sh - submit ONE algorithm for BI
# Usage: bash hpc/cross-dataset-evaluation/scripts/submit_bi.sh mcsplit

ALGO=$1
if [ -z "$ALGO" ]; then
    echo "Usage: bash hpc/cross-dataset-evaluation/scripts/submit_bi.sh <algo>"
    echo "Algos: mcsplit rl ll dal dsb rrsplit symsplit"
    exit 1
fi

mkdir -p hpc/logs hpc/cross-dataset-evaluation/scripts hpc/cross-dataset-evaluation/results

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
SCRIPT="hpc/cross-dataset-evaluation/scripts/job_${ALGO}_bi.slurm"
sed "s/{ALGO}/$ALGO/g; s/{NJOBS}/$NJOBS/g; s/{DATASET}/bi/g" \
    hpc/cross-dataset-evaluation/scripts/job_array.slurm > "$SCRIPT"
sed -i 's/\r$//' "$SCRIPT"

# Verify replacement worked
if grep -q "{DATASET}" "$SCRIPT"; then
    echo "ERROR: {DATASET} was not replaced in $SCRIPT"
    exit 1
fi

echo "Submitting $ALGO BI: $((NJOBS+1)) jobs"
sbatch "$SCRIPT"