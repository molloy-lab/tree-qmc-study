#!/bin/bash
#SBATCH --qos=highmem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64gb
#SBATCH --account=cbcb
#SBATCH --partition=cbcb
#SBATCH --exclusive
#SBATCH --time=24:00:00
#SBATCH --constraint=EPYC-7313

lscpu
cat /proc/meminfo

# Define directories and input
GROUPDIR="/fs/cbcb-lab/ekmolloy"
PROJDIR="$GROUPDIR/ekmolloy/tree-qmc-study/han2024wtreeqmc"
DATADIR="$PROJDIR/data/morel2022asteroid-plants"
GTRE_FILE="all-relabeled-raxml-ng.LG+G.geneTree.newick"
MYINFO="plants"

echo "Submitting TREE-QMC job..."
srun \
    --nodes 1 \
    --ntasks 1 \
    --cpus-per-task 32 \
    --mem=36gb \
    ./b5_run_treeqmc.sh $MYINFO $DATADIR $GTRE_FILE  
wait
echo "...done"

echo "Submitting Asteroid job..."
srun \
    --nodes 1 \
    --ntasks 1 \
    --cpus-per-task 32 \
    --mem=36gb \
    ./c5_run_asteroid.sh $MYINFO $DATADIR $GTRE_FILE 
wait
echo "...done"

echo "Submitting Asteroid job..."
srun \
    --nodes 1 \
    --ntasks 1 \
    --cpus-per-task 32 \
    --mem=36gb \
    ./d5_run_wastrid.sh $MYINFO $DATADIR $GTRE_FILE 
wait
echo "...done"


echo "Submitting ASTER job..."
srun \
    --nodes 1 \
    --ntasks 1 \
    --cpus-per-task 32 \
    --mem=36gb \
    ./e5_run_aster.sh $MYINFO $DATADIR $GTRE_FILE 
wait
echo "...done"


scontrol show job $SLURM_JOB_ID
scontrol write batch_script $SLURM_JOB_ID -

