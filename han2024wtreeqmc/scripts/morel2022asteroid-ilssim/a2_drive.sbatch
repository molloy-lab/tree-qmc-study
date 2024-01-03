#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16gb
#SBATCH --partition=cbcb
#SBATCH --account=cbcb
#SBATCH --qos=highmem
#SBATCH --time=00:04:00

MYINFO="SIM,MISS,S,F,SITES,MODL,BL,D,L,T,GC,P,POP,MS,MF,SEED"
MYMODL=$(echo $MODL | sed 's/_/,/g')  # Line in CSV

# Define directories
GROUPDIR="/fs/cbcb-lab/ekmolloy"
PROJDIR="$GROUPDIR/ekmolloy/tree-qmc-study/han2024wtreeqmc"

# Define true species tree
INDIR="$GROUPDIR/group/data/morel2022asteroid/simulations/withils/$MODL"
STRE_TRUE="$INDIR/species_trees/speciesTree.newick"

# Define gene trees
OUTDIR="$PROJDIR/data/morel2022asteroid-ilssim/$MODL"
mkdir -p $OUTDIR

cp $STRE_TRUE $OUTDIR

GTRE_TRUE="$OUTDIR/true_gene_trees.tre"
GTRE_ESTI="$OUTDIR/raxmlng_gene_trees.tre"
if [ ! -e $GTRE_TRUE ]; then
    cat $INDIR/families/family_*/gene_trees/true.true.geneTree.newick | sed 's/;/;\n/g' | sed 's/_0_0//g' > $GTRE_TRUE
fi
if [ ! -e $GTRE_ESTI ]; then
    cat $INDIR/families/family_*/gene_trees/raxml-ng.GTR+G.geneTree.newick | sed 's/_0_0//g' > $GTRE_ESTI
fi

./a2_get_properties.sh $MYMODL $OUTDIR $STRE_TRUE $GTRE_TRUE $GTRE_ESTI

