#!/bin/bash

# https://slurm.schedmd.com/srun.html
# https://rc-docs.northeastern.edu/en/latest/using-discovery/srun.html

#SBATCH --qos=throughput
#SBATCH --exclude=heron[00-11]
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=36gb
#SBATCH --exclusive
#SBATCH --time=04:00:00

lscpu
cat /proc/meminfo

echo "Submitting job..."
srun \
    --nodes 1 \
    --ntasks 1 \
    --cpus-per-task 12 \
    --mem=36gb \
    ./g3_compute_quartet_score.sh $SCAL $NGEN $NBPS $REPL
wait
echo "...done"

scontrol show job $SLURM_JOB_ID
scontrol write batch_script $SLURM_JOB_ID -

