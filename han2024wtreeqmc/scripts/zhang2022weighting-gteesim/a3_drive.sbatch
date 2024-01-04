#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16gb
#SBATCH --partition=cbcb
#SBATCH --account=cbcb
#SBATCH --qos=highmem
#SBATCH --time=00:04:00

MYMODL="$REPL,$NBPS"

# Define directories
GROUPDIR="/fs/cbcb-lab/ekmolloy"
PROJDIR="$GROUPDIR/ekmolloy/tree-qmc-study/han2024wtreeqmc"

# Define true species tree
INDIR1="$GROUPDIR/group/data/zhang2022weighting-dryad/S100/$REPL"
INDIR2="$GROUPDIR/group/data/zhang2018astral3/ASTRALIII/S100/$REPL"
STRE_TRUE="$INDIR1/s_tree.trees"
GTRE_TRUE="$INDIR2/truegenetrees"
GTRE_ESTI="$INDIR1/bestMLestimatedgenetree/estimatedgenetre_${NBPS}.gtr.rerooted.final.contracted.non"

# Define gene trees
OUTDIR="$PROJDIR/data/zhang2022weighting-gteesim/S100/$REPL"

./a3_get_properties.sh $MYMODL $OUTDIR $STRE_TRUE $GTRE_TRUE $GTRE_ESTI $NBPS
