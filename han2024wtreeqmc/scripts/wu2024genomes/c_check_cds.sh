#!/bin/bash

exit

GROUPDIR="/fs/cbcb-lab/ekmolloy"
PROJDIR="/fs/cbcb-lab/ekmolloy/ekmolloy/tree-qmc-study/han2024wtreeqmc"
GENEDIR="$GROUPDIR/group/data/wu2024genomes/1_besttrees_cds5756"
OUTDIR="$PROJDIR/data/wu2024genomes/gene-trees-abayes"

cd $GENEDIR

mkdir 1-check-align_OG0003316_cds
mv *align_OG0003316_cds* 1-check-align_OG0003316_cds

mkdir 2-check-align_OG0006453_cds
mv *align_OG0006453_cds* 2-check-align_OG0006453_cds

mkdir 3-check-align_OG0009507_cds
mv *align_OG0009507_cds* 3-check-align_OG0009507_cds

mkdir 4-check-align_OG0009765_cds
mv *align_OG0009765_cds* 4-check-align_OG0009765_cds

mkdir 5-check-align_OG0009997_cds
mv *align_OG0009997_cds* 5-check-align_OG0009997_cds

cat *check*/RAxML_bestTree* > check_cds_raxml.tre
mv check_cds_raxml.tre $OUTDIR/..

cd $OUTDIR

mkdir 1-check-align_OG0003316_cds
mv *align_OG0003316_cds* 1-check-align_OG0003316_cds

mkdir 2-check-align_OG0006453_cds
mv *align_OG0006453_cds* 2-check-align_OG0006453_cds

mkdir 3-check-align_OG0009507_cds
mv *align_OG0009507_cds* 3-check-align_OG0009507_cds

mkdir 4-check-align_OG0009765_cds
mv *align_OG0009765_cds* 4-check-align_OG0009765_cds

mkdir 5-check-align_OG0009997_cds
mv *align_OG0009997_cds* 5-check-align_OG0009997_cds

cat *check*/*treefile > check_cds_raxml_abayes.tre
mv check_cds_raxml_abayes.tre ..

