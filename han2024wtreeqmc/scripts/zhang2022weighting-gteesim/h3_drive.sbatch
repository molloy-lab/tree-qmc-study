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
        "aster_h_t16" \
        "wastrid_s" \
        "wtreeqmc_wf_n2" \
        "wtreeqmc_wh_n0" \
        "wtreeqmc_wh_n1" \
        "wtreeqmc_wh_n2" \
        "wtreeqmc_wn_n2" \
        "wtreeqmc_ws_n2" )

echo "Running branch_support estimation..."
for MTHD in ${MTHDS[@]}; do
    bash h3_compute_quartet_score.sh  $REPL $NBPS $SUPP $NGEN $MTHD
done
echo "done"

