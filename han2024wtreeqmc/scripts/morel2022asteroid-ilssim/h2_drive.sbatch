#!/bin/bash
#SBATCH --qos=highmem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64gb
#SBATCH --account=cbcb
#SBATCH --partition=cbcb
#SBATCH --time=00:30:00
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

echo "Running branch support analyses..."
bash h2_compare_branch_support.sh $MYMODL $OUTDIR
echo "done"

