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

MYINFO="SIM,MISS,S,F,SITES,MODL,BL,D,L,T,GC,P,POP,MS,MF,SEED"
MYMODL=$(echo $MODL | sed 's/_/,/g')  # Line in CSV

# Define directories
GROUPDIR="/fs/cbcb-lab/ekmolloy"
PROJDIR="$GROUPDIR/ekmolloy/tree-qmc-study/han2024wtreeqmc"

# Define true species tree
INDIR="$GROUPDIR/group/data/morel2022asteroid/simulations/withils/$MODL"
STRE_TRUE="$INDIR/species_trees/speciesTree.newick"

# Define input gene trees
OUTDIR="$PROJDIR/data/morel2022asteroid-ilssim/$MODL"
GTRE_FILE="raxmlng_gene_trees.tre"

echo "Running wASTRID..."
bash b2_run_wastrid.sh $MYMODL $OUTDIR $GTRE_FILE $STRE_TRUE
echo "...done"

echo "Running Asteroid..."
bash c2_run_asteroid.sh $MYMODL $OUTDIR $GTRE_FILE $STRE_TRUE
echo "...done"

echo "Running wTREE-QMC..."
bash d2_run_wtreeqmc.sh $MYMODL $OUTDIR $GTRE_FILE $STRE_TRUE 
echo "...done"

echo "Running ASTRAL-III"
bash e2_run_astral3.sh $MYMODL $OUTDIR $GTRE_FILE $STRE_TRUE
echo "...done"

echo "Running ASTER"
bash f2_run_aster.sh $MYMODL $OUTDIR $GTRE_FILE $STRE_TRUE
echo "done"

