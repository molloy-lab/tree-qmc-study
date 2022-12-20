#!/bin/bash

# https://slurm.schedmd.com/srun.html
# https://rc-docs.northeastern.edu/en/latest/using-discovery/srun.html

#SBATCH --qos=throughput
#SBATCH --exclude=cbcbgpu00,cbcbgpu02,ibis[02-06,08-09,11-18,20-25],locust,redbud,spruce,tern[00-01],gum
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
    ./b1_run_treeqmc.sh $NTAX $PSIZE $SRATE $REPL $GTRE $NGEN
wait
echo "...done"

scontrol show job $SLURM_JOB_ID
scontrol write batch_script $SLURM_JOB_ID -

