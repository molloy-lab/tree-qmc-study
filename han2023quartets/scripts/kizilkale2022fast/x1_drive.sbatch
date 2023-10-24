#!/bin/bash

# https://slurm.schedmd.com/srun.html
# https://rc-docs.northeastern.edu/en/latest/using-discovery/srun.html

#SBATCH --qos=throughput
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16gb
#SBATCH --time=01:00:00

lscpu
cat /proc/meminfo

MYINFO="SIMNO,S,M,H,MINVAF,ISAV,N,FP,FN,NA,D,L"
MYMODL=$(echo $MODL | sed 's/_//g' | sed 's/-/,/g')  # Line in CSV

# Define directories
GROUPDIR="/fs/cbcb-lab/ekmolloy"
PROJECTDIR="$GROUPDIR/ekmolloy/species-to-tumor-study"
DATDIR="$PROJECTDIR/data/kizilkale2022fast/$MODL"


echo "Mapping mutations to trees..."
srun \
    --nodes 1 \
    --ntasks 1 \
    --cpus-per-task 1 \
    --mem=16gb \
    ./x1_root_tree.sh $MYMODL $DATDIR
wait
echo "...done"


scontrol show job $SLURM_JOB_ID
scontrol write batch_script $SLURM_JOB_ID -

