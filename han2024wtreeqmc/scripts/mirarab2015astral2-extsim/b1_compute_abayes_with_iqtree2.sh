#!/bin/bash

NTAX=$1  # Number of taxa
HGHT=$2  # Species tree height (number of generations)
RATE=$3  # Speciation rate
REPL=$4  # Replicate number
GNUM=$5  # Gene number

MODL="model.$NTAX.$HGHT.$RATE"

# Define directories
LABDIR="/fs/cbcb-lab/ekmolloy"
PROJDIR="$LABDIR/ekmolloy/tree-qmc-study/han2024wtreeqmc"

# Define software files
IQTREE2="$LABDIR/group/software/iqtree-2.2.2.6-Linux/bin/iqtree2"
GETNEWICK="$PROJDIR/tools/get_newick_from_stacked_file.py"
GETPHYLIP="$PROJDIR/tools/get_phylip_from_stacked_file.py"

# Define input files
INDIR="$LABDIR/group/data/mirarab2015astral2/$MODL/$REPL"
ALGN="all-genes.phylip"
GTRE="estimatedgenetre"

if [ ! -e $INDIR/$ALGN ]; then
    echo "$INDIR/$ALGN does not exist!"
    exit 1
fi

if [ ! -e $INDIR/$GTRE ]; then
    echo "$INDIR/$GTRE does not exist!"
    exit 1
fi

# Define output files
OUTDIR="$PROJDIR/data/mirarab2015astral2-extsim/$MODL/$REPL"
GTAG=$(seq -f "%04g" $GNUM $GNUM )

cd $OUTDIR

if [ -e $GTRE.abayes-$GTAG ]; then
    echo "$OUTDIR/$GTRE.abayes-$GTAG already exists!"
    exit 1
fi

# 1 - split the phylip file so each gene has an alignment
python3 $GETPHYLIP \
    -i $INDIR/$ALGN \
    -n $GNUM \
    -o tmp-$ALGN-$GTAG

# 2 - split gene tree file so each gene has a tree
python3 $GETNEWICK \
    -i $INDIR/$GTRE \
    -n $GNUM \
    -o tmp-$GTRE-$GTAG

# 3 - estimate abayes support for each gene tree with iqtree2
$IQTREE2 \
    -s tmp-$ALGN-$GTAG \
    -te tmp-$GTRE-$GTAG \
    -m GTR+I+G4 \
    -abayes \
    -pre tmp-$GTRE.abayes-$GTAG \
    -T 1

sed 's/\///g' tmp-$GTRE.abayes-$GTAG.treefile > $GTRE.abayes-$GTAG

# 4 - clean up
rm tmp-$ALGN-$GTAG
rm tmp-$GTRE-$GTAG
rm tmp-$GTRE.abayes-$GTAG.*

