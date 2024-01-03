#!/bin/bash

MODL="ssim_veryhighmiss_s25_f1000_sites100_GTR_bl1.0_d0.0_l0.0_t0.0_gc0.0_p0.0_pop50000000_ms0.6_mf0.6_seed3000"
MODL="ssim_veryhighmiss_s125_f1000_sites100_GTR_bl1.0_d0.0_l0.0_t0.0_gc0.0_p0.0_pop50000000_ms0.6_mf0.6_seed3036"
echo "Processing $MODL"

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

#echo "Running wASTRID..."
#bash b2_run_wastrid.sh $MYMODL $OUTDIR $GTRE_FILE $STRE_TRUE
#echo "...done"

#echo "Running Asteroid..."
#bash c2_run_asteroid.sh $MYMODL $OUTDIR $GTRE_FILE $STRE_TRUE
#echo "...done"

echo "Running wTREE-QMC..."
bash d2_run_wtreeqmc.sh $MYMODL $OUTDIR $GTRE_FILE $STRE_TRUE 
echo "...done"

#echo "Running ASTRAL-III"
#bash e2_run_astral3.sh $MYMODL $OUTDIR $GTRE_FILE $STRE_TRUE
#echo "...done"

#echo "Running ASTER"
#bash f2_run_aster.sh $MYMODL $OUTDIR $GTRE_FILE $STRE_TRUE
#echo "done"

