#!/bin/bash
#SBATCH --job-name=mcs_{ALGO}
#SBATCH --array=0-{NJOBS}%500
#SBATCH --time=00:12:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=1
#SBATCH --output=hpc/logs/{ALGO}_%A_%a.out
#SBATCH --error=hpc/logs/{ALGO}_%A_%a.err

ALGO={ALGO}
INSTANCE_LIST=hpc/instances_{ALGO}.txt
INSTANCE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" $INSTANCE_LIST)
A=$(echo $INSTANCE | awk '{print $1}')
B=$(echo $INSTANCE | awk '{print $2}')

bash hpc/run_instance.sh $ALGO $A $B