#!/bin/bash

# https://slurm.schedmd.com/srun.html
# https://rc-docs.northeastern.edu/en/latest/using-discovery/srun.html

#SBATCH --qos=throughput
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=36gb
#SBATCH --exclusive
#SBATCH --time=18:00:00

lscpu
cat /proc/meminfo

echo "Submitting job..."
srun \
    --nodes 1 \
    --ntasks 1 \
    --cpus-per-task 32 \
    --mem=36gb \
    ./e2_run_extract_and_weight_quartets.sh $SCAL $NGEN $NBPS $REPL
wait
echo "...done"

scontrol show job $SLURM_JOB_ID
scontrol write batch_script $SLURM_JOB_ID -

