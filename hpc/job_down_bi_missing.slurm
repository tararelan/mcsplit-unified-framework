#!/bin/bash
#SBATCH --job-name=mcs_down
#SBATCH --partition=large-short
#SBATCH --array=0-11%168
#SBATCH --time=00:40:00
#SBATCH --mem=40G
#SBATCH --cpus-per-task=1
#SBATCH --output=hpc/logs/down_%A_%a.out
#SBATCH --error=hpc/logs/down_%A_%a.err

ALGO=down
DATASET=bi
WITH_REF=1
SWAP=0
INSTANCE_LIST=hpc/missing/down_bi.txt
INSTANCE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" $INSTANCE_LIST)
A=$(echo $INSTANCE | awk '{print $1}')
B=$(echo $INSTANCE | awk '{print $2}')

bash hpc/run_instance.sh $ALGO $A $B $DATASET $WITH_REF $SWAP