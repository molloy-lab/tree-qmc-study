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

MTHDS=( "true_stree" \
        "aster_h" \
        "wtreeqmc_wh_n2" )

echo "Running branch_support estimation..."
for MTHD in ${MTHDS[@]}; do
    for SUPP in "abayes" "bs"; do
        bash h3_compute_quartet_score.sh  $REPL $NBPS $SUPP $NGEN $MTHD
    done
done
echo "done"

