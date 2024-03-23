import copy
import pandas
import numpy
import sys

"""
https://stackoverflow.com/questions/556405/what-do-real-user-and-sys-mean-in-the-output-of-time1
Real is wall clock time - time from start to finish of the call. This is all elapsed time including time slices used by other processes and time the process spends blocked (for example if it is waiting for I/O to complete).
User is the amount of CPU time spent in user-mode code (outside the kernel) within the process. This is only actual CPU time used in executing the process. Other processes and time the process spends blocked do not count towards this figure.
Sys is the amount of CPU time spent in the kernel within the process. This means executing CPU time spent in system calls within the kernel, as opposed to library code, which is still running in user-space. Like 'user', this is only CPU time used by the process. See below for a brief description of kernel mode (also known as 'supervisor' mode) and the system call mechanism.
"""

#sys.exit("DONE RUNNING")

def reformat_timing(data):
    [m, s] = data.split('m')
    m = float(m)
    s = float(s.replace('s', ''))

    totals = m * 60.0 + s
    totalm = totals / 60.0
    totalh = totalm / 60.0

    return [totals, totalm, totalh]

ste_df = pandas.read_csv("all_species_tree_error.csv", keep_default_na=False)
mrt_df = pandas.read_csv("all_runtime.csv", keep_default_na=False)

cols = ["NTAX", "NGEN", "NBPS", "BLSC", "PSIZ", "MISS",
        "REPL", "MTHD",
        "SEFN", "SEFNR",
        "SEFP", "SEFPR",
        "NODE", "SECS"]

mthds = ["asteroid",
         "aster_v1.16.3.4",
         "wastrid_vanilla",
         "wtreeqmc_wf_n0",
         "wtreeqmc_wf_n0_refined",
         "wtreeqmc_wf_n1",
         "wtreeqmc_wf_n1_refined",
         "wtreeqmc_wf_n1_shared",
         "wtreeqmc_wf_n1_shared_refined",
         "wtreeqmc_wf_n2",
         "wtreeqmc_wf_n2_refined",
         "wtreeqmc_wf_n2_shared",
         "wtreeqmc_wf_n2_shared_refined",
         "wtreeqmc_wn_n2",
         "wtreeqmc_wn_n2_refined"]

index = 1
experiments = ["varypsiz", "varyntax", "varyngen",
               "varynbps", "varyblsc", "varymiss"]
for do in experiments:
    rows = []

    ss = ["s50"]
    fs = ["f1000"]
    sitess = ["sites100"]
    bls = ["bl1.0"]
    pops = ["pop50000000"]
    mss = ["ms0.6"]  # always same as mf
    repls = ["seed30" + str(x).zfill(2) for x in range(0, 50)]
    
    if do == "varypsiz":
        #for varying pop size 5 * 50 = 250
        pops = ["pop10", "pop50000000", "pop100000000", "pop500000000", "pop1000000000"]
    elif do == "varyntax":
        # for varying number of taxa, 5 * 50 = 250
        ss = ["s25", "s75", "s50", "s100", "s125", "s150"]  # dup s50
    elif do == "varyngen":
        # for varying number of genes, 3 * 50 = 150
        fs = ["f250", "f500", "f1000", "f2000"]  # dup f1000
    elif do == "varynbps":
        # for varying number of sites 3 * 50 = 150
        sitess = ["sites50", "sites100", "sites200", "sites500"]  # dup sites100
    elif do == "varyblsc":
        # for varying branch length scaler 5 * 50 = 250
        bls = ["bl0.05", "bl0.1", "bl1.0", "bl10.0", "bl100.0", "bl200.0"]  # dup "bl1.0"
    elif do == "varymiss":
        # for varying something 5 * 50 = 250
        mss = ["ms0.5", "ms0.55", "ms0.6", "ms0.65", "ms0.7", "ms0.75"]  # dup "ms0.6"
    else:
        sys.exit()

    rows = []

    for s in ss:
        for f in fs:
            for sites in sitess:
                for bl in bls:
                    for pop in pops:
                        for ms in mss:
                            for repl in repls:
                                print("%d %s %s %s %s %s %s %s %s" % (index, do, s, f, sites, bl, pop, ms, repl))

                                for mthd in mthds:
                                    xste_df = ste_df[(ste_df["S"] == s) &
                                                     (ste_df["F"] == f) &
                                                     (ste_df["SITES"] == sites) &
                                                     (ste_df["BL"] == bl) &
                                                     (ste_df["POP"] == pop) &
                                                     (ste_df["MS"] == ms) &
                                                     (ste_df["SEED"] == repl) &
                                                     (ste_df["MTHD"] == mthd)]

                                    if xste_df.shape[0] != 1:
                                        sys.exit("  1 ERROR - %s!\n" % mthd)

                                    if (xste_df.FN.values[0] == '') or (xste_df.FN.values[0] == "NA"):
                                        sefn = "NA"
                                        sefp = "NA"
                                        sefnr = "NA"
                                        sefpr = "NA"
                                    else:
                                        sefn = int(xste_df.FN.values[0])
                                        sefp = int(xste_df.FP.values[0])
                                        sefnr = float(sefn) / float(xste_df.I1.values[0])
                                        sefpr = float(sefp) / float(xste_df.I2.values[0])

                                    xmrt_df = mrt_df[(mrt_df["S"] == s) &
                                                     (mrt_df["F"] == f) &
                                                     (mrt_df["SITES"] == sites) &
                                                     (mrt_df["BL"] == bl) &
                                                     (mrt_df["POP"] == pop) &
                                                     (mrt_df["MS"] == ms) &
                                                     (mrt_df["SEED"] == repl) &
                                                     (mrt_df["MTHD"] == mthd.replace("_refined", ''))]

                                    if xmrt_df.shape[0] != 1:
                                        sys.exit("  2 RUNTIME - %s!\n" % mthd)

                                    if (xmrt_df.real.values[0] == '') or (xmrt_df.real.values[0] == "NA"):
                                        secs = "NA"
                                        node = "NA"
                                    else:
                                        mrt = reformat_timing(xmrt_df.real.values[0])
                                        secs = mrt[0]

                                        node = xmrt_df.NODE.values[0]
                                        if node.find("cbcb") < 0:
                                            sys.exit("Method not run on CBCB node - check architecture!\n")

                                    row = {}
                                    row["NTAX"] = int(s.replace('s', ''))
                                    row["NGEN"] = int(f.replace('f', ''))
                                    row["NBPS"] = int(sites.replace('sites', ''))
                                    row["BLSC"] = float(bl.replace('bl', ''))
                                    row["PSIZ"] = int(pop.replace('pop', ''))
                                    row["MISS"] = float(ms.replace('ms', ''))
                                    row["REPL"] = int(repl.replace('seed30', '')) + 1
                                    row["MTHD"] = mthd
                                    row["NODE"] = node
                                    row["SEFN"] = sefn
                                    row["SEFNR"] = sefnr
                                    row["SEFP"] = sefp
                                    row["SEFPR"] = sefpr
                                    row["SECS"] = secs
                                    rows.append(row)

                                index += 1

    df = pandas.DataFrame(rows, columns=cols)
    df.to_csv("data-" + do + "-error-and-timings.csv",
              sep=',', na_rep="NA",header=True, index=False)
