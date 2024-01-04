#!/bin/bash

REPL=$1
NBPS=$2
SUPP=$3
NGEN=$4

#echo "Running wASTRID..."
#bash b3_run_wastrid.sh $REPL $NBPS $SUPP $NGEN
#echo "...done"

echo "Running Asteroid..."
bash c3_run_asteroid.sh $REPL $NBPS $SUPP $NGEN
echo "...done"

echo "Running wTREE-QMC-h..."
bash d3_run_wtreeqmc.sh $REPL $NBPS $SUPP $NGEN
echo "...done"

echo "Running TREE-QMC again..."
bash e3_run_wtreeqmc_wn.sh $REPL $NBPS $SUPP $NGEN
echo "done"

echo "Running TREE-QMC again again..."
bash f3_run_wtreeqmc_varyn.sh $REPL $NBPS $SUPP $NGEN
echo "done"

#echo "Running ASTER-h..."
#bash g3_run_asterh.sh $REPL $NBPS $SUPP $NGEN
#echo "...done"


