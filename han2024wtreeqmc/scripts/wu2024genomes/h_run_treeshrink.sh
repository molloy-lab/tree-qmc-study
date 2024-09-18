#!/bin/bash

exit

GROUPDIR="/fs/cbcb-lab/ekmolloy"
PROJDIR="/fs/cbcb-lab/ekmolloy/ekmolloy/tree-qmc-study/han2024wtreeqmc"
DATADIR="$PROJDIR/data/wu2024genomes/gene-trees-abayes"

cd $DATADIR

run_treeshrink.py -t ./cds_raxml_abayes.tre -m "per-species" -q 0.05
mv cds_raxml_abayes_treeshrink/output.tre cds_raxml_abayes_treeshrink.tre
mv cds_raxml_abayes_treeshrink/output.txt cds_raxml_abayes_treeshrink.txt
mv cds_raxml_abayes_treeshrink/output_summary.txt cds_raxml_abayes_treeshrink_summary.txt
rm -rf cds_raxml_abayes_treeshrink

run_treeshrink.py -t ./introns_raxml_abayes.tre -m "per-species" -q 0.05
mv introns_raxml_abayes_treeshrink/output.tre introns_raxml_abayes_treeshrink.tre
mv introns_raxml_abayes_treeshrink/output.txt introns_raxml_abayes_treeshrink.txt
mv introns_raxml_abayes_treeshrink/output_summary.txt introns_raxml_abayes_treeshrink_summary.txt
rm -rf introns_raxml_abayes_treeshrink

cat cds_raxml_abayes.tre \
    introns_raxml_abayes.tre \
    > cds_and_introns_raxml_abayes.tre

cat cds_raxml_abayes_treeshrink.tre \
    introns_raxml_abayes_treeshrink.tre \
    > cds_and_introns_raxml_abayes_treeshrink.tre

