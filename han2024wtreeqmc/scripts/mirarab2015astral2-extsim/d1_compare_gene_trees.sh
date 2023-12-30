#!/bin/bash

# IMPORTANT OBSERVATION:
# Running IQTREE refines polytomies from FastTree (for identical sequences) and 
# does NOT assign support to those refined branches

NTAX=$1  # Number of taxa
HGHT=$2  # Species tree height (number of generations)
RATE=$3  # Speciation rate
REPL=$4  # Replicate number

MODL="model.$NTAX.$HGHT.$RATE"

# Define directories
LABDIR="/fs/cbcb-lab/ekmolloy"
PROJDIR="$LABDIR/ekmolloy/tree-qmc-study/han2024wtreeqmc"

# Define software files
CMPTREES="$PROJDIR/tools/compare_two_tree_lists.py"
EXBRANCH="$PROJDIR/tools/extract_gene_tree_branch_info.py"

# Define input files
OLDDIR="$LABDIR/group/data/mirarab2015astral2/$MODL/$REPL"
NEWDIR="$PROJDIR/data/mirarab2015astral2-extsim/$MODL/$REPL"

TRUE_GTRE="$OLDDIR/truegenetrees"
ESTI_GTRE_SH="$OLDDIR/estimatedgenetre"
ESTI_GTRE_AB="$NEWDIR/estimatedgenetre.abayes"

if [ ! -e $TRUE_GTRE ]; then
    echo "$TRUE_GTRE does not exist!"
    exit 1
fi

if [ ! -e $ESTI_GTRE_SH ]; then
    echo "$ESTI_GTRE_SH does not exist!"
    exit 1
fi

if [ ! -e $ESTI_GTRE_AB ]; then
    echo "$ESTI_GTRE_AB does not exist!"
    exit 1
fi

CSVF="$NEWDIR/true_vs_sh_gene_trees.csv"
if [ ! -e $CSVF ]; then
    python3 $EXBRANCH \
        -t $TRUE_GTRE \
        -e $ESTI_GTRE_SH \
        -p "$NTAX,$HGHT,$RATE,$REPL,SH" \
        &> $CSVF
fi

CSVF="$NEWDIR/true_vs_abayes_gene_trees.csv"
if [ ! -e $CSVF ]; then
    python3 $EXBRANCH \
        -t $TRUE_GTRE \
        -e $ESTI_GTRE_AB \
        -p "$NTAX,$HGHT,$RATE,$REPL,AB" \
        -o $ESTI_GTRE_AB.fixed \
        &> $CSVF
fi

CSVF="$NEWDIR/sh_vs_abayes_gene_trees.csv"
if [ ! -e $CSVF ]; then
    python3 $CMPTREES \
        -l1 $ESTI_GTRE_SH \
        -l2 $ESTI_GTRE_AB.fixed \
	-p "$NTAX,$HGHT,$RATE,$REPL" \
       &> $CSVF	
fi

if [ -e $CSVF ]; then
    MAXFN=$(sed 's/,/ /g' $CSVF | awk '{print $11}' | sort -n | tail -n1)
    if [ $MAXFN -ne 0 ]; then
        echo "Bad gene tree!" > $NEWDIR/bad_gene_tree.txt
    else
        rm $NEWDIR/estimatedgenetre.abayes-*
    fi
fi

