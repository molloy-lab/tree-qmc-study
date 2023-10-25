import pandas
import numpy
import sys

inputs = ["all_cell_lineage_tree_error.csv",
          "all_mut_pair_error.csv",
          "all_runtime.csv",
          "all_cell_lineage_tree_quartet_score.csv"]

for inputf in inputs:
    base = ["SIMNO", "S", "M", "H", "MINVAF", "ISAV", "N", "ALPHA_FP", "BETA_FN", "GAMMA_NA", "D", "L"]
    if inputf == "all_cell_lineage_tree_error.csv":
        outputf = "all_cell_lineage_tree_error_fixed.csv"
        header = base + ["MTHD", "NL", "I1", "I2", "FN", "FP", "RF"]
    elif inputf == "all_mut_pair_error.csv":
        outputf = "all_mut_pair_error_fixed.csv"
        header = base + ["MTHD", 
                         "SL_TP", "SL_TN", "SL_FP", "SL_FN",
                         "DL_TP", "DL_TN", "DL_FP", "DL_FN",
                         "AD_TP", "AD_TN", "AD_FP", "AD_FN", "AD_FLIP"]
    elif inputf == "all_runtime.csv":
        outputf = "all_runtime_fixed.csv"
        header = base + ["MTHD", "NODE", "NTHREADS", "real", "user", "sys"]
    else:
        outputf = "all_cell_lineage_tree_quartet_score_fixed.csv"
        header = base + ["MTHD", "QS", "NQS"]

    ncol2 = len(header)
    ncol1 = ncol2 - 6

    with open (inputf, 'r') as fin, open(outputf, 'w') as fout:
        fout.write(','.join(header) + '\n')
        for line in fin:
            words = line.split(',')
            nword = len(words)
            if nword == ncol1:
                simno = words[0]
                n = words[1]
                m = words[2]
                fp = words[3]
                fn = words[4]
                na = words[5]
                end = ','.join(words[6:])
                new_line = str("%s,NA,%s,NA,NA,NA,%s,%s,%s,%s,NA,NA,%s\n" % (simno, m, n, fp, fn, na, end))
                fout.write(new_line)
            elif nword == ncol2:
                fout.write(line)
            else:
                sys.exit("Unrecognized number of elements: %s" % line)

