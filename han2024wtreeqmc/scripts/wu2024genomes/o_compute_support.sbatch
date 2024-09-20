#!/bin/bash
#SBATCH --qos=highmem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64gb
#SBATCH --account=cbcb
#SBATCH --partition=cbcb
#SBATCH --exclusive
#SBATCH --time=72:00:00
#SBATCH --constraint=EPYC-7313

lscpu
cat /proc/meminfo

# Define directories
GROUPDIR="/fs/cbcb-lab/ekmolloy"
PROJDIR="/fs/cbcb-lab/ekmolloy/ekmolloy/tree-qmc-study/han2024wtreeqmc"
GENEDIR="$PROJDIR/data/wu2024genomes/gene-trees-abayes"
OUTDIR="$PROJDIR/data/wu2024genomes/species-trees"
ASTERH="$GROUPDIR/group/software-compiled-on-EPYC-7313/ASTER/bin/astral-hybrid"
# --mode	Integer	1	1: hybrid weighting, 2: support only, 3: length only, 4: unweighted
#--scoring	Preset		Scoring the full species tree file after `-c` without exploring other topologies (`-r 1 -s 0`)

cd $OUTDIR

# TREE A
STRE="rooted_treeA_wastrid_support_none_nonefilter"
GTRE="cds_and_introns_raxml_abayes"
$ASTERH --mode 4 \
        -t 16 -u 2 -r 1 -s 0 \
        --scoring --root allgtrmss \
        -c $STRE.tre \
        -i $GENEDIR/${GTRE}.tre \
        -o ${STRE}_w_unweighted_support.tre \
        &> ${STRE}_w_unweighted_support.log

STRE="rooted_treeA_wastrid_support_none_nonefilter"
GTRE="cds_and_introns_raxml_abayes_treeshrink"
$ASTERH --mode 4 \
        -t 16 -u 2 -r 1 -s 0 \
        --scoring --root allgtrmss \
        -c $STRE.tre \
        -i $GENEDIR/${GTRE}.tre \
        -o ${STRE}_w_unweighted_support_after_treeshrink.tre \
        &> ${STRE}_w_unweighted_support_after_treeshrink.log

# TREE B
STRE="rooted_treeB_wtreeqmc_hybrid_length"
GTRE="cds_and_introns_raxml_abayes"
$ASTERH --mode 1 --bayes \
        -t 16 -u 2 -r 1 -s 0 \
        --scoring --root allgtrmss \
        -c $STRE.tre \
        -i $GENEDIR/${GTRE}.tre \
        -o ${STRE}_w_hybrid_support.tre \
        &> ${STRE}_w_hybrid_support.log

# TREE C
STRE="rooted_treeC_wtreeqmc_support_none_nonefilter"
GTRE="cds_and_introns_raxml_abayes"
$ASTERH --mode 4 \
        -t 16 -u 2 -r 1 -s 0 \
        --scoring --root allgtrmss \
        -c $STRE.tre \
        -i $GENEDIR/${GTRE}.tre \
        -o ${STRE}_w_unweighted_support.tre \
        &> ${STRE}_w_unweighted_support.log

STRE="rooted_treeC_wtreeqmc_support_none_nonefilter"
GTRE="cds_and_introns_raxml_abayes_treeshrink"
$ASTERH --mode 4 \
        -t 16 -u 2 -r 1 -s 0 \
        --scoring --root allgtrmss \
        -c $STRE.tre \
        -i $GENEDIR/${GTRE}.tre \
        -o ${STRE}_w_unweighted_support_after_treeshrink.tre \
        &> ${STRE}_w_unweighted_support_after_treeshrink.log

# TREE D
STRE="rooted_treeD_waster_hybrid"
GTRE="cds_and_introns_raxml_abayes"
$ASTERH --mode 1 --bayes \
        -t 16 -u 2 -r 1 -s 0 \
        --scoring --root allgtrmss \
        -c $STRE.tre \
        -i $GENEDIR/${GTRE}.tre \
        -o ${STRE}_w_hybrid_support.tre \
        &> ${STRE}_w_hybrid_support.log

# TREE E
STRE="rooted_treeE_waster_length"
GTRE="cds_and_introns_raxml_abayes"
$ASTERH --mode 3 \
        -t 16 -u 2 -r 1 -s 0 \
        --scoring --root allgtrmss \
        -c $STRE.tre \
        -i $GENEDIR/${GTRE}.tre \
        -o ${STRE}_w_length_support.tre \
        &> ${STRE}_w_length_support.log

# TREE F
STRE="rooted_treeF_waster_support_none"
GTRE="cds_and_introns_raxml_abayes"
$ASTERH --mode 4 \
        -t 16 -u 2 -r 1 -s 0 \
        --scoring --root allgtrmss \
        -c $STRE.tre \
        -i $GENEDIR/${GTRE}.tre \
        -o ${STRE}_w_unweighted_support.tre \
        &> ${STRE}_w_unweighted_support.log

STRE="rooted_treeF_waster_support_none"
GTRE="cds_and_introns_raxml_abayes_treeshrink"
$ASTERH --mode 4 \
        -t 16 -u 2 -r 1 -s 0 \
        --scoring --root allgtrmss \
        -c $STRE.tre \
        -i $GENEDIR/${GTRE}.tre \
        -o ${STRE}_w_unweighted_support_after_treeshrink.tre \
        &> ${STRE}_w_unweighted_support_after_treeshrink.log

# TREE G
STRE="rooted_treeG_waster_nonefilter"
GTRE="cds_and_introns_raxml_abayes_treeshrink"
$ASTERH --mode 4 \
        -t 16 -u 2 -r 1 -s 0 \
        --scoring --root allgtrmss \
        -c $STRE.tre \
        -i $GENEDIR/${GTRE}.tre \
        -o ${STRE}_w_unweighted_support_after_treeshrink.tre \
        &> ${STRE}_w_unweighted_support_after_treeshrink.log

#rooted_treeA_wastrid_support_none_nonefilter_w_unweighted_support_after_treeshrink.log:Score: 5.159284703e+10
#rooted_treeA_wastrid_support_none_nonefilter_w_unweighted_support.log:Score: 5.621946598e+10
#rooted_treeB_wtreeqmc_hybrid_length_w_hybrid_support.log:Score: 3.621396671e+10
#rooted_treeC_wtreeqmc_support_none_nonefilter_w_unweighted_support_after_treeshrink.log:Score: 5.162360162e+10
#rooted_treeC_wtreeqmc_support_none_nonefilter_w_unweighted_support.log:Score: 5.625060623e+10
#rooted_treeD_waster_hybrid_w_hybrid_support.log:Score: 3.624623623e+10
#rooted_treeE_waster_length_w_length_support.log:Score: 3.785523365e+10
#rooted_treeF_waster_support_none_w_unweighted_support_after_treeshrink.log:Score: 5.166484467e+10
#rooted_treeF_waster_support_none_w_unweighted_support.log:Score: 5.629419268e+10
#rooted_treeG_waster_nonefilter_w_unweighted_support_after_treeshrink.log:Score: 5.166663719e+10
