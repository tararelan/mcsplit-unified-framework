#!/bin/bash
#SBATCH --job-name=mcs_dsb
#SBATCH --partition=large-short
#SBATCH --array=0-2679%168
#SBATCH --time=00:40:00
#SBATCH --mem=40G
#SBATCH --cpus-per-task=1
#SBATCH --output=hpc/logs/dsb_%A_%a.out
#SBATCH --error=hpc/logs/dsb_%A_%a.err

ALGO=dsb
DATASET=bi
WITH_REF=1
SWAP=0
INSTANCE_LIST=hpc/missing/dsb_bi.txt
INSTANCE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" $INSTANCE_LIST)
A=$(echo $INSTANCE | awk '{print $1}')
B=$(echo $INSTANCE | awk '{print $2}')

bash hpc/run_instance.sh $ALGO $A $B $DATASET $WITH_REF $SWAP