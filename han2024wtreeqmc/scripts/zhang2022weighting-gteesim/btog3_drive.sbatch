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

echo "Running wASTRID..."
bash b3_run_wastrid.sh $REPL $NBPS $SUPP $NGEN
echo "...done"

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

echo "Running ASTER-h with 1 threads..."
bash g3_run_asterh.sh $REPL $NBPS $SUPP $NGEN
echo "...done"

