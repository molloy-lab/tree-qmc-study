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

echo "Running branch support estimation..."
MTHDS=( "true_stree" "aster_h_t16" "wtreeqmc_wh_n2" )
for MTHD in ${MTHDS[@]}; do
    bash i1_compute_quartet_score.sh $NTAX $HGHT $RATE $REPL $SUPP $NGEN $MTHD
done
echo "done"

