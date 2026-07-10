#!/bin/bash
# submit_bi.sh — submit ONE algorithm for BI
# Usage: bash submit_bi.sh mcsplit

ALGO=$1
if [ -z "$ALGO" ]; then
    echo "Usage: bash submit_bi.sh <algo>"
    echo "Algos: mcsplit rl dal ll dsb rrsplit symsplit"
    exit 1
fi

mkdir -p hpc/logs hpc/cross-dataset-evaluation

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
SCRIPT="hpc/cross-dataset-evaluation/job_${ALGO}_bi.slurm"
sed "s/{ALGO}/$ALGO/g; s/{NJOBS}/$NJOBS/g" hpc/cross-dataset-evaluation/job_array.slurm > "$SCRIPT"
sed -i 's/\r$//' "$SCRIPT"

echo "Submitting $ALGO BI: $((NJOBS+1)) jobs"
sbatch "$SCRIPT"