#!/bin/bash
#SBATCH --qos=highmem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64gb
#SBATCH --account=cbcb
#SBATCH --partition=cbcb
#SBATCH --time=24:00:00
#SBATCH --constraint=EPYC-7313

lscpu
cat /proc/meminfo

echo "Running TREE-QMC again..."
bash i1_run_wtreeqmc_wn.sh  $NTAX $HGHT $RATE $REPL $SUPP $NGEN
echo "done"

echo "Running TREE-QMC again again..."
bash j1_run_wtreeqmc_varyn.sh  $NTAX $HGHT $RATE $REPL $SUPP $NGEN
echo "done"

