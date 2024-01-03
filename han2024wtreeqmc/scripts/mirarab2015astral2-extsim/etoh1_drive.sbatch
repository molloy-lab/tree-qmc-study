#!/bin/bash
#SBATCH --qos=highmem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64gb
#SBATCH --account=cbcb
#SBATCH --partition=cbcb
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

echo "Running ASTER-h with 16 threads..."
bash h1_run_asterh_t16.sh  $NTAX $HGHT $RATE $REPL $SUPP $NGEN
echo "...done"

