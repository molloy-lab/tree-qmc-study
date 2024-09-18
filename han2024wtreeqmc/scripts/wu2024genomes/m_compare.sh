#!/bin/bash

GROUPDIR="/fs/cbcb-lab/ekmolloy"
PROJDIR="/fs/cbcb-lab/ekmolloy/ekmolloy/tree-qmc-study/han2024wtreeqmc"
OUTDIR="$PROJDIR/data/wu2024genomes/species-trees"
COMPARE="$PROJDIR/tools/compare_two_trees.py"

cd $OUTDIR


# ASTRID
TREE1="wastrid_length_cds_and_introns_raxml_abayes.tre"
TREE2="wastrid_support_cds_and_introns_raxml_abayes.tre"
TREE3="wastrid_none_cds_and_introns_raxml_abayes.tre"
TREE4="wastrid_none_cds_and_introns_raxml_abayes_treeshrink.tre"

python3 $COMPARE -t1 $TREE1 -t2 $TREE2  # length vs. support
# 125,125,125,122,122,18,18

python3 $COMPARE -t1 $TREE2 -t2 $TREE3  # support vs none
#125,125,125,122,122,0,0

python3 $COMPARE -t1 $TREE2 -t2 $TREE4  # support vs none w/filter taxa
#125,125,125,122,122,0,0

TREEA=$TREE2 # ASTRID, everything but length <- plot (b)
TREEX=$TREE1 # ASTRID length <- don't plot


# TREE-QMC
TREE5="wtreeqmc_hybrid_n2_cds_and_introns_raxml_abayes.tre"
TREE6="wtreeqmc_length_n2_cds_and_introns_raxml_abayes.tre"
TREE7="wtreeqmc_support_n2_cds_and_introns_raxml_abayes.tre"
TREE8="wtreeqmc_none_n2_cds_and_introns_raxml_abayes.tre"
TREE9="wtreeqmc_none_n2_cds_and_introns_raxml_abayes_treeshrink.tre"

python3 $COMPARE -t1 $TREE5 -t2 $TREE6 # hybrid vs. length
# 125,125,125,122,122,0,0

python3 $COMPARE -t1 $TREE5 -t2 $TREE7 # hybrid vs. support
# 125,125,125,122,122,3,3

python3 $COMPARE -t1 $TREE7 -t2 $TREE8 # support vs. none
# 125,125,125,122,122,0,0

python3 $COMPARE -t1 $TREE7 -t2 $TREE9 # support vs. none w/filter
# 125,125,125,122,122,0,0


TREEB=$TREE5  # TREEQMC hybrid, length <- plot (c)
TREEC=$TREE7  # TREEQMC support, none, none w/ filter <- plot (d)

python3 $COMPARE -t1 $TREEA -t2 $TREEB # ASTRID vs. TREE-QMC hybrid/length
#125,125,125,122,122,3,3

python3 $COMPARE -t1 $TREEA -t2 $TREEC # ASTRID vs. TREE-QMC support, none, none w/filter
#125,125,125,122,122,4,4


# ASTRAL
TREE10="waster_hybrid_t16_cds_and_introns_raxml_abayes.tre"
TREE11="waster_length_t16_cds_and_introns_raxml_abayes.tre"
TREE12="waster_support_t16_cds_and_introns_raxml_abayes.tre"
TREE13="waster_none_t16_cds_and_introns_raxml_abayes.tre"
TREE14="waster_none_t16_cds_and_introns_raxml_abayes_treeshrink.tre"

python3 $COMPARE -t1 $TREE12 -t2 $TREE13  # support vs. none
#125,125,125,122,122,0,0

python3 $COMPARE -t1 $TREE12 -t2 $TREE14  # support/none vs. none w/filter <- more sensitive to long branches?
#125,125,125,122,122,4,4

python3 $COMPARE -t1 $TREE10 -t2 $TREE12  # hybrid vs. support/none
#125,125,125,122,122,6,6

python3 $COMPARE -t1 $TREE10 -t2 $TREE14  # hybrid vs. none w/filter taxa
#125,125,125,122,122,3,3

python3 $COMPARE -t1 $TREE10 -t2 $TREE11  # hybrid vs. length
#125,125,125,122,122,4,4

python3 $COMPARE -t1 $TREE11 -t2 $TREE14  # length vs. none w/filter taxa <- length is similar to filtering
#125,125,125,122,122,1,1

TREED=$TREE10 # ASTRAL hybrid <- plot (c)
TREEE=$TREE11 # ASTRAL length
TREEF=$TREE12 # ASTRAL support
TREEG=$TREE14 # ASTRAL none w/filter <- plot (e)

python3 $COMPARE -t1 $TREEA -t2 $TREED # ASTRID vs. ASTRAL hybrid
# 125,125,125,122,122,14,14

python3 $COMPARE -t1 $TREEX -t2 $TREED # ASTRID length vs. ASTRAL hybrid
# 125,125,125,122,122,22,22

python3 $COMPARE -t1 $TREEA -t2 $TREEB # ASTRID vs. TREE-QMC hybrid/length
#125,125,125,122,122,3,3

python3 $COMPARE -t1 $TREEB -t2 $TREED # TREE-QMC hybrid/length vs. ASTRAL hybrid
# 125,125,125,122,122,14,14

python3 $COMPARE -t1 $TREEB -t2 $TREEE # TREE-QMC hybrid/length vs. ASTRAL length
# 125,125,125,122,122,12,12

python3 $COMPARE -t1 $TREEB -t2 $TREEF # TREE-QMC hybrid/length vs. ASTRAL  support/none
# 125,125,125,122,122,14,14

python3 $COMPARE -t1 $TREEB -t2 $TREEG # TREE-QMC hybrid/length vs. ASTRAL none w/filter
#125,125,125,122,122,13,13


# ASTEROID
TREE15="asteroid_cds_and_introns_raxml_abayes.tre.bestTree.newick"
TREE16="asteroid_cds_and_introns_raxml_abayes_treeshrink.tre.bestTree.newick"

python3 $COMPARE -t1 $TREE15 -t2 $TREE16
#125,125,125,122,122,4,4

TREEH=$TREE15 # Asteroid
TREEI=$TREE16 # Asteroid after treeshrink

python3 $COMPARE -t1 $TREEA -t2 $TREEH # ASTRID vs. Asteroid 
#125,125,125,122,122,4,4

python3 $COMPARE -t1 $TREEA -t2 $TREEI # ASTRID vs. Asteroid w/ filter
#125,125,125,122,122,0,0

