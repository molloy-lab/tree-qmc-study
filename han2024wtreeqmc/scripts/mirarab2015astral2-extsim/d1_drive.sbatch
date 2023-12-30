#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16gb
#SBATCH --partition=cbcb
#SBATCH --account=cbcb
#SBATCH --qos=highmem
#SBATCH --time=24:00:00

./d1_compare_gene_trees.sh $NTAX $HGHT $RATE $REPL

