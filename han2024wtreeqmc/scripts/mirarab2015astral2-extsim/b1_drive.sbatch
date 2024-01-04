#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8gb
#SBATCH --partition=cbcb
#SBATCH --account=cbcb
#SBATCH --qos=highmem
#SBATCH --time=24:00:00

for GNUM in `seq $GS $GE`; do
    echo "Processing gene $GNUM..."
    ./b1_compute_abayes_with_iqtree2.sh $NTAX $HGHT $RATE $REPL $GNUM
    echo "...done"
done

