There are three collection of data sets, referred to as 
+ `mirarab2015astral2`, 
+ `mahbub2021wqfm-aviansim`, and 
+ `mahbub2021wqfm-mammaliansim`. 

All of these data sets contain true and estimated species trees, true and estimated gene trees, and refined gene trees (i.e., gene trees with polytomies arbitrarily refinedusing TREE-QMC). These trees were used to compute different quantities (e.g. error), which are in the CSVS files. The trees have been uploaded to dryad and the CSV files are available here.

# Trees for mirarab2015astral2
Tree files are organized with the following structure `<MODL>/<REPL>`, where 
+ `MODL` denotes the model parameters and 
+ `REPL` denotes the replicate data set number (1-50). 

The model parameters have the form `model.<NTAX>.<STRHT>.<SRATE>`, where 
+ `NTAX` is the number of taxa (`10`, `50`, `100`, `500`, `1000`)
+ `STRHT` is the species tree height (5X = `10000000`, 1X = `2000000`, 0.5X = `500000`)
+ `SRATE` is the speciation rate (deep = `0.0000001` and shallow = `0.000001`)

Note that true species trees, true gene trees, and estimated gene trees were downloaded from [here](https://sites.google.com/eng.ucsd.edu/datasets/astral/astral-ii) in September 2021.

# Trees for mahbub2021wqfm-aviansim
Tree files are organized with the following structure `<MODL>/<REPL>`, where 
+ `MODL` denotes the model parameters and 
+ `REPL` denotes the replicate data set number (R1-R20). 

The model parameters have the form `<SCAL>-<NGEN>-<NBPS>`, where 
+ `SCAL` is the species tree scale (`0.5X`, `1X`, `2X`)
+ `NGEN` is the number of gene trees
+ `NBPS` is indicates whether gene trees are `true` or estimated (sequence length: `500`)

The data files indicate whether the gene trees are true or estmated.

Note that 

# Trees for mahbub2021wqfm-mammaliansim
Tree files are organized with the following structure `<MODL>/<REPL>`, where 
+ `MODL` denotes the model parameters and 
+ `REPL` denotes the replicate data set number (R1-R20). 

The model parameters have the form `<SCAL>.<NGEN>.<NBPS>`, where 
+ `SCAL` is the species tree scale (0.5X = `scale2d`, 1X = `noscale`, 2X = `scale2u`)
+ `NGEN` is the number of gene trees
+ `NBPS` is indicates whether gene trees are `true` or estimated (sequence length: `250`, `500`, `1000`, `1500`)

# CSV for species tree estimation
Three CSV files contain information about species tree estimation and thus include the method `MTHD` used and the input it was given (note: CSVS files include the model parameters indicated on each row). 

The header of **all_quartet_score.csv** ends in `GTRE,NGEN,MTHD,QS,NQS`
+ `QS` is the quartet score where
+ `NQS` is the normalized quartet score

The header of **all_quartets_runtime.csv** and **all_runtime.csv** ends in `NODE,real,user,sys` where
+ `NODE` is teh compute node used to run the analysis
+ `real` is the wallclock time
+ `user` and `sys` are not used in the analyses

The header of **all_species_tree_error.csv** ends in `NL,I1,I2,FN,FP,RF` where
+ `NL` is the number of leaves
+ `I1` is the number of internal branches in tree 1 (true species tree)
+ `I2` is the number of internal branches in tree 2 (species tree estimated with MTHD)
+ `FN` is the number of false negative branches (i.e. number of branches in true 1 that are missing from tree 2)
+ `FP` is the number of false positive branches (i.e. number of branches in tree 2 that are missing from tree 1)
+ `RF` is the normalized Robinson-Foulds distance between tree 1 and tree 2 `FN + FP / (2*NL - 6)`

# CSV files for data set properties
Three CSV files contain information about the trees and thus the indicate the `GENE` number (which corresponds to line that gene tree is on in the file).

The header of **true_vs_estimated_gene_trees.csv** ends in `NL,I1,I2,FN,FP,RF`.
The quantities are similar to *all_species_tree_error.csv* except that tree 1 is the true gene tree (on line `GENE`) and tree 2 is the estimated gene tree (on line `GENE`). Thus, this file contains information about amount of gene tree estimation error (GTEE).

The header of **true_species_tree_vs_true_gene_trees.csv** ends in `NL,I1,I2,FN,FP,RF`. Again, these quantiteis are similar to *all_species_tree_error.csv* except that tree 1 is the true species tree and tree 2 is the true gene tree (on line `GENE`). Thus, this file contains information about the amount of incomplete lineage sorting (ILS).

The header of **true_species_tree_vs_estimated_gene_trees.csv** ends in `NL,I1,I2,FN,FP,RF`. Again, these quantiteis are similar to **all_species_tree_error.csv** except that tree 1 is the true species tree and tree 2 is the estimated gene tree (on line `GENE`). Thus, this files contains information about of gene tree heterogeneity due to either GTEE or ILS.
