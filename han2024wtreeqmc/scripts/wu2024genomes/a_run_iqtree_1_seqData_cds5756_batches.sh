#!/bin/bash

GROUPDIR="/fs/cbcb-lab/ekmolloy"
PROJDIR="/fs/cbcb-lab/ekmolloy/ekmolloy/tree-qmc-study/han2024wtreeqmc"

IQTREE="$GROUPDIR/software/iqtree-2.3.5-Linux-intel/bin/iqtree2"
OUTDIR="$PROJDIR/data/wu2024genomes/gene-trees-abayes/1_seqData_cds5756"
DATADIR="$GROUPDIR/group/data/wu2024genomes/1_seqData_cds5756"
GENEDIR="$GROUPDIR/group/data/wu2024genomes/1_besttrees_cds5756"

source "a_run_iqtree_1_seqData_cds5756.sh"

