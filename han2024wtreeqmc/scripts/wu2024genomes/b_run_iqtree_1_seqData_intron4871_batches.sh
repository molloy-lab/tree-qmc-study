#!/bin/bash

GROUPDIR="/fs/cbcb-lab/ekmolloy"
PROJDIR="/fs/cbcb-lab/ekmolloy/ekmolloy/tree-qmc-study/han2024wtreeqmc"

IQTREE="$GROUPDIR/software/iqtree-2.3.5-Linux-intel/bin/iqtree2"
OUTDIR="$PROJDIR/data/wu2024genomes/gene-trees-abayes/1_besttrees_intron4871"
DATADIR="$GROUPDIR/group/data/wu2024genomes/1_seqData_intron4871"
GENEDIR="$GROUPDIR/group/data/wu2024genomes/1_besttrees_intron4871"

source "b_run_iqtree_1_seqData_intron4871.sh"

