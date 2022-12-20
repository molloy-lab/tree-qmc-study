#!/bin/bash

# https://slurm.schedmd.com/srun.html
# https://rc-docs.northeastern.edu/en/latest/using-discovery/srun.html

#SBATCH --qos=throughput
#SBATCH --exclude=cbcbgpu00,cbcbgpu02,ibis[02-06,08-09,11-18,20-25],locust,redbud,spruce,tern[00-01],gum
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16gb
#SBATCH --exclusive
#SBATCH --time=12:00:00

srun \
    --nodes 1 \
    --ntasks 1 \
    --cpus-per-task 1 \
    --mem=16gb \
    ./a1_get_properties.sh $NTAX $PSIZE $SRATE
wait

