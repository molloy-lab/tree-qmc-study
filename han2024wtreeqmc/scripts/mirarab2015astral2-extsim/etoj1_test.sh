#!/bin/bash

NTAX=$1
HGHT=$2
RATE=$3
REPL=$4
SUPP=$5
NGEN=$6

echo "Running wASTRID..."
bash e1_run_wastrid.sh  $NTAX $HGHT $RATE $REPL $SUPP $NGEN
echo "...done"

echo "Running Asteroid..."
bash f1_run_asteroid.sh  $NTAX $HGHT $RATE $REPL $SUPP $NGEN
echo "...done"

echo "Running wTREE-QMC-h..."
bash g1_run_wtreeqmc.sh  $NTAX $HGHT $RATE $REPL $SUPP $NGEN
echo "...done"

echo "Running ASTER-h..."
bash h1_run_asterh.sh  $NTAX $HGHT $RATE $REPL $SUPP $NGEN
echo "...done"

echo "Running TREE-QMC again..."
bash i1_run_wtreeqmc_again.sh  $NTAX $HGHT $RATE $REPL $SUPP $NGEN
echo "done"

echo "Running TREE-QMC again again..."
bash j1_run_wtreeqmc_again_again.sh  $NTAX $HGHT $RATE $REPL $SUPP $NGEN
echo "done"

