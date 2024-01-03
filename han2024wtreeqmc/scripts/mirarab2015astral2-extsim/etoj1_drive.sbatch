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

