import pandas
import numpy
import sys

"""
https://stackoverflow.com/questions/556405/what-do-real-user-and-sys-mean-in-the-output-of-time1
Real is wall clock time - time from start to finish of the call. This is all elapsed time including time slices used by other processes and time the process spends blocked (for example if it is waiting for I/O to complete).
User is the amount of CPU time spent in user-mode code (outside the kernel) within the process. This is only actual CPU time used in executing the process. Other processes and time the process spends blocked do not count towards this figure.
Sys is the amount of CPU time spent in the kernel within the process. This means executing CPU time spent in system calls within the kernel, as opposed to library code, which is still running in user-space. Like 'user', this is only CPU time used by the process. See below for a brief description of kernel mode (also known as 'supervisor' mode) and the system call mechanism.
"""

def reformat_timing(data):
    [m, s] = data.split('m')
    m = float(m)
    s = float(s.replace('s', ''))

    totals = m * 60.0 + s
    totalm = totals / 60.0
    totalh = totalm / 60.0

    return [totals, totalm, totalh]


ste_df = pandas.read_csv("all_cell_lineage_tree_error_fixed.csv", keep_default_na=False)
mpe_df = pandas.read_csv("all_mut_pair_error_fixed.csv", keep_default_na=False)
mrt_df = pandas.read_csv("all_runtime_fixed.csv", keep_default_na=False)
nqs_df = pandas.read_csv("all_cell_lineage_tree_quartet_score_fixed.csv", keep_default_na=False)

cols = ["NCELL", "NMUT", "ALPHA_FP","BETA_FN","GAMMA_NA",
        "S", "H", "MINVAF", "ISAV", "D", "L", "REPL",
        "MTHD", "NODE", "SECS",
        "NLEAF", "NINT_ESTI", "SEFP", "NINT_TRUE", "SEFN",
        "SLP", "SLR", "DLP", "DLR", "ADP", "ADR", "ADFLIP",
        "QS", "NQS"]

rows = []

mthds1 = ["huntress_v0.1.2.0_default",
          "scistree_v1.2.0.6",
          "fastral_wroot",
          "fastme_v2.1.5_wroot",
          "treeqmcbip_v1.0.0_n2_wroot"]

mthds2 = ["huntress_v0.1.2.0_default",
          "scistree_v1.2.0.6_wmuts",
          "fastral_wrootx_wmuts",
          "fastme_v2.1.5_wrootx_wmuts",
          "treeqmcbip_v1.0.0_n2_wrootx_wmuts"]

mthds3 = ["huntress_v0.1.2.0_default",
          "scistree_v1.2.0.6",
          "fastral_wrootx",
          "fastme_v2.1.5_wrootx",
          "treeqmcbip_v1.0.0_n2_wrootx"]

# Extended figure 4 (matches 5 and 6) -  80 / 470 model conditions
ncxms = [["n300", "m300"], ["n300", "m1000"], ["n1000", "m300"]]

fns = ["fn0.05", "fn0.2"]                         # False negative rate
fns = ["fn0.2"]                                   # False negative rate
fp = "fp0.001"                                    # No false positives
na = "na0.05"                                     # No missing data
s = "s100"                                        # 
h = "h1"                                          #
minvaf = "minVAF0.005"                            #
isav = "ISAV0"                                    # ISA violated false
d = "d0"                                          # No doublets
l = "l1000000"                                    #
repls = ["simNo" + str(x) for x in range(1, 11)]  # Different trees (so diff number's of internal branches)

index = 0

for ncxm in ncxms:
    ncell = ncxm[0]
    nmut = ncxm[1]
    for fn in fns:
        xste_df = ste_df[(ste_df["N"] == ncell) &
                         (ste_df["M"] == nmut) &
                         (ste_df["BETA_FN"] == fn) &
                         (ste_df["ALPHA_FP"] == fp) &
                         (ste_df["GAMMA_NA"] == na) &
                         (ste_df["S"] == s) &
                         (ste_df["H"] == h) &
                         (ste_df["MINVAF"] == minvaf) &
                         (ste_df["ISAV"] == isav) &
                         (ste_df["D"] == d) &
                         (ste_df["L"] == l) &
                         (ste_df["MTHD"] == "scistree_v1.2.0.6")]

        tmp1 = xste_df.I1.values
        tmp2 = xste_df.I2.values

        print(numpy.mean(tmp1))
        print(numpy.mean(tmp2))

        for repl in repls:
            print("%d %s %s %s %s" % (index, ncell, nmut, fn, repl))

            for i in range(len(mthds1)):
                mthd1 = mthds1[i]
                mthd2 = mthds2[i]
                mthd3 = mthds3[i]

                # Process cell lineage tree error
                xste_df = ste_df[(ste_df["N"] == ncell) &
                                 (ste_df["M"] == nmut) &
                                 (ste_df["BETA_FN"] == fn) &
                                 (ste_df["ALPHA_FP"] == fp) &
                                 (ste_df["GAMMA_NA"] == na) &
                                 (ste_df["S"] == s) &
                                 (ste_df["H"] == h) &
                                 (ste_df["MINVAF"] == minvaf) &
                                 (ste_df["ISAV"] == isav) &
                                 (ste_df["D"] == d) &
                                 (ste_df["L"] == l) &
                                 (ste_df["SIMNO"] == repl) &
                                 (ste_df["MTHD"] == mthd2)]

                if xste_df.shape[0] != 1:
                    sys.exit("  2 ERROR - %s!\n" % mthd2)

                if (xste_df.FN.values[0] == '') or (xste_df.FN.values[0] == "NA"):
                    nl = "NA"
                    ni1 = "NA"
                    ni2 = "NA"
                    sefn = "NA"
                    sefp = "NA"
                else:
                    nl = int(xste_df.NL.values[0])
                    ni1 = int(xste_df.I1.values[0])
                    ni2 = int(xste_df.I2.values[0])
                    sefn = int(xste_df.FN.values[0])
                    sefp = int(xste_df.FP.values[0])

                # Process mut pair error
                xmpe_df = mpe_df[(mpe_df["N"] == ncell) &
                                 (mpe_df["M"] == nmut) &
                                 (mpe_df["BETA_FN"] == fn) &
                                 (mpe_df["ALPHA_FP"] == fp) &
                                 (mpe_df["GAMMA_NA"] == na) &
                                 (mpe_df["S"] == s) &
                                 (mpe_df["H"] == h) &
                                 (mpe_df["MINVAF"] == minvaf) &
                                 (mpe_df["ISAV"] == isav) &
                                 (mpe_df["D"] == d) &
                                 (mpe_df["L"] == l) &
                                 (mpe_df["SIMNO"] == repl) &
                                 (mpe_df["MTHD"] == mthd2)]

                if xmpe_df.shape[0] != 1:
                    print(xmpe_df)
                    sys.exit("  2 ERROR - %s!\n" % mthd2)
                    
                if (xmpe_df.SL_TP.values[0] == '') or (xmpe_df.SL_TP.values[0] == "NA"):
                    sl_p = "NA"
                    sl_r = "NA"
                    dl_p = "NA"
                    dl_r = "NA"
                    ad_p = "NA"
                    ad_r = "NA"
                    ad_flip_rate = "NA"
                else:
                    sl_tp = int(xmpe_df.SL_TP.values[0])
                    sl_fp = int(xmpe_df.SL_FP.values[0])
                    sl_fn = int(xmpe_df.SL_FN.values[0])
                    sl_p = sl_tp / (sl_tp + sl_fp)
                    sl_r = sl_tp / (sl_tp + sl_fn)

                    dl_tp = int(xmpe_df.DL_TP.values[0])
                    dl_fp = int(xmpe_df.DL_FP.values[0])
                    dl_fn = int(xmpe_df.DL_FN.values[0])
                    dl_p = dl_tp / (dl_tp + dl_fp)
                    dl_r = dl_tp / (dl_tp + dl_fn)

                    ad_tp = int(xmpe_df.AD_TP.values[0])
                    ad_fp = int(xmpe_df.AD_FP.values[0])
                    ad_fn = int(xmpe_df.AD_FN.values[0])
                    ad_flip = int(xmpe_df.AD_FLIP.values[0])
                    ad_p = ad_tp / (ad_tp + ad_fp + ad_flip)
                    ad_r = ad_tp / (ad_tp + ad_fn + ad_flip)
                    ad_flip_rate = ad_flip / (ad_tp + ad_fp + ad_flip)

                # Process runtime
                xmrt_df = mrt_df[(mrt_df["N"] == ncell) &
                                 (mrt_df["M"] == nmut) &
                                 (mrt_df["BETA_FN"] == fn) &
                                 (mrt_df["ALPHA_FP"] == fp) &
                                 (mrt_df["GAMMA_NA"] == na) &
                                 (mrt_df["S"] == s) &
                                 (mrt_df["H"] == h) &
                                 (mrt_df["MINVAF"] == minvaf) &
                                 (mrt_df["ISAV"] == isav) &
                                 (mrt_df["D"] == d) &
                                 (mrt_df["L"] == l) &
                                 (mrt_df["SIMNO"] == repl) &
                                 (mrt_df["MTHD"] == mthd1)]

                if xmrt_df.shape[0] != 1:
                    sys.exit("  3 ERROR - %s!\n" % mthd1)

                if (xmrt_df.real.values[0] == '') or (xmrt_df.real.values[0] == "NA"):
                    secs = "NA"
                    node = "NA"
                else:
                    mrt = reformat_timing(xmrt_df.real.values[0])
                    secs = mrt[0]
                    node = xmrt_df.NODE.values[0]
                    if (mthd1 == "huntress_v0.1.2.0_default"):
                        savenode = node
                    else:
                        if node != savenode:
                            sys.exit("Methods run on different nodes!\n")


                # Process quartet score
                xnqs_df = nqs_df[(nqs_df["N"] == ncell) &
                                 (nqs_df["M"] == nmut) &
                                 (nqs_df["BETA_FN"] == fn) &
                                 (nqs_df["ALPHA_FP"] == fp) &
                                 (nqs_df["GAMMA_NA"] == na) &
                                 (nqs_df["S"] == s) &
                                 (nqs_df["H"] == h) &
                                 (nqs_df["MINVAF"] == minvaf) &
                                 (nqs_df["ISAV"] == isav) &
                                 (nqs_df["D"] == d) &
                                 (nqs_df["L"] == l) &
                                 (nqs_df["SIMNO"] == repl) &
                                 (nqs_df["MTHD"] == mthd3)]

                if xnqs_df.shape[0] != 1:
                    sys.exit("  4 QUARTET SCORE - %s!\n" % mthd3)

                if (xnqs_df.QS.values[0] == '') or (xnqs_df.QS.values[0] == "NA"):
                    qscore = "NA"
                    nqscore = "NA"
                else:
                    qscore = int(xnqs_df.QS.values[0])
                    nqscore = float(xnqs_df.NQS.values[0])

                # Build data frame
                row = {}
                row["NCELL"] = ncell
                row["NMUT"] = nmut
                row["ALPHA_FP"] = fp
                row["BETA_FN"] = fn
                row["GAMMA_NA"] = na
                row["S"] = s
                row["H"] = h
                row["MINVAF"] = minvaf
                row["ISAV"] = isav
                row["D"] = d
                row["L"] = l
                row["REPL"] = repl
                row["MTHD"] = mthd2
                row["NODE"] = node
                row["SECS"] = secs
                row["NLEAF"] = nl
                row["NINT_TRUE"] = ni1
                row["NINT_ESTI"] = ni2
                row["SEFN"] = sefn
                row["SEFP"] = sefp
                row["SLP"] = sl_p
                row["SLR"] = sl_r
                row["DLP"] = dl_p
                row["DLR"] = dl_r
                row["ADP"] = ad_p
                row["ADR"] = ad_r
                row["ADFLIP"] = ad_flip_rate
                row["QS"] = qscore
                row["NQS"] = nqscore
                rows.append(row)

                index += 1

df = pandas.DataFrame(rows, columns=cols)
df.to_csv("data-all-error-and-timings_exfig4-6.csv",
          sep=',', na_rep="NA",header=True, index=False)
