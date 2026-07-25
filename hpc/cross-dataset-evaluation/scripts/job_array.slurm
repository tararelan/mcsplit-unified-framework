#!/bin/bash
#SBATCH --job-name=mcs_{ALGO}_{DATASET}
#SBATCH --partition=large-short
#SBATCH --array=0-{NJOBS}%168
#SBATCH --time=00:40:00
#SBATCH --mem=36G
#SBATCH --cpus-per-task=1
#SBATCH --output=hpc/logs/{ALGO}_{DATASET}_%A_%a.out
#SBATCH --error=hpc/logs/{ALGO}_{DATASET}_%A_%a.err

ALGO={ALGO}
DATASET={DATASET}
INSTANCE_LIST=hpc/instances_{ALGO}_{DATASET}.txt
INSTANCE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" $INSTANCE_LIST)
A=$(echo $INSTANCE | awk '{print $1}')
B=$(echo $INSTANCE | awk '{print $2}')

# bash hpc/cross-dataset-evaluation/scripts/run_instance.sh $ALGO $A $B $DATASET

bash hpc/cross-dataset-evaluation/scripts/run_instance_reference.sh $ALGO $A $B $DATASET
