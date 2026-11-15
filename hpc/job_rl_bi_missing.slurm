#!/bin/bash
#SBATCH --job-name=mcs_rl
#SBATCH --partition=large-short
#SBATCH --array=0-69%168
#SBATCH --time=00:40:00
#SBATCH --mem=40G
#SBATCH --cpus-per-task=1
#SBATCH --output=hpc/logs/rl_%A_%a.out
#SBATCH --error=hpc/logs/rl_%A_%a.err

ALGO=rl
DATASET=bi
WITH_REF=1
SWAP=0
INSTANCE_LIST=hpc/missing/rl_bi.txt
INSTANCE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" $INSTANCE_LIST)
A=$(echo $INSTANCE | awk '{print $1}')
B=$(echo $INSTANCE | awk '{print $2}')

bash hpc/run_instance.sh $ALGO $A $B $DATASET $WITH_REF $SWAP