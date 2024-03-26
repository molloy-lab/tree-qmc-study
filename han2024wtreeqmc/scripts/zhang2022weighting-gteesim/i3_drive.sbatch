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

NBPSS=( 200 400 800 1600 )
NGENS=( 50 200 500 1000 )
SUPPS=( "abayes" "bs" )
for NBPS in ${NBPSS[@]}; do
    for SUPP in ${SUPPS[@]}; do
        for NGEN in ${NGENS[@]}; do
            echo "Working on $NBPS bps $NGEN gen"
            bash i3_compare_branch_support.sh  $REPL $NBPS $SUPP $NGEN
        done
    done
done

