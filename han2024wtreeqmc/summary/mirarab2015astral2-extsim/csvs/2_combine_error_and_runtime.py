import csv
import pandas
import numpy
import sys

"""
https://stackoverflow.com/questions/556405/what-do-real-user-and-sys-mean-in-the-output-of-time1
Real is wall clock time - time from start to finish of the call. This is all elapsed time including time slices used by other processes and time the process spends blocked (for example if it is waiting for I/O to complete).
User is the amount of CPU time spent in user-mode code (outside the kernel) within the process. This is only actual CPU time used in executing the process. Other processes and time the process spends blocked do not count towards this figure.
Sys is the amount of CPU time spent in the kernel within the process. This means executing CPU time spent in system calls within the kernel, as opposed to library code, which is still running in user-space. Like 'user', this is only CPU time used by the process. See below for a brief description of kernel mode (also known as 'supervisor' mode) and the system call mechanism.
"""

"""
Previous studies removed:
+ remove 1 replicate (i.e. 41) from 10-taxon
+ remove 3 replicates (i.e. 21 and 41) from 50-taxon (next is 27 with 590ish genes)
+ remove 2 replicates (i.e. 8 and 47) from 100-taxon
+ remove 3 replicates (i.e. 8, 15, 49) from 200 taxon/500 K/1e-06

Failures due to wall clock time for single-threaded ASTER-h:
+ remove 4 replicates (6, 8, 33, and 38) from 1000-taxon
"""

sys.exit("DONE RUNNING")

namemap = {}
namemap["caml"] = "CA-ML"
namemap["treeqmc_n2_v1.0.0"] = "TQMC-n0-origstudy"
namemap["aster_h_t16"] = "ASTER-wh (16 threads)"
namemap["aster_h"] = "ASTER-wh (1 thread)"
namemap["wastrid_s"] = "ASTRID-ws"
namemap["wtreeqmc_wn_n2"] = "TQMC-n2"
namemap["wtreeqmc_ws_n2"] = "TQMC-ws_n2"
namemap["wtreeqmc_wh_n2"] = "TQMC-wh_n2"
namemap["wtreeqmc_wh_n1"] = "TQMC-wh_n1"
namemap["wtreeqmc_wh_n0"] = "TQMC-wh_n0"

def reformat_timing(data):
    [m, s] = data.split('m')
    m = float(m)
    s = float(s.replace('s', ''))

    totals = m * 60.0 + s
    totalm = totals / 60.0
    totalh = totalm / 60.0

    return [totals, totalm, totalh]


df1 = pandas.read_csv("all_species_tree_error.csv", keep_default_na=False)
df2 = pandas.read_csv("aster_h_t1_species_tree_error.csv", keep_default_na=False)
ste_df = pandas.concat([df1, df2])

df1 = pandas.read_csv("all_runtime.csv", keep_default_na=False)
df2 = pandas.read_csv("aster_h_t1_runtime.csv", keep_default_na=False)
mrt_df = pandas.concat([df1, df2])

cols = ["NTAX", "ILSL", "SPEC", "REPL", "SUPP", "NGEN", 
        "MTHD", "NODE", "SEFN", "SEFNR", "SECS"]

ngens = [1000]
supps = ["sh"]

for do in ["varyils", "varyntax"]:
    if do == "varyntax":
        sys.stdout.write("Increasing number of taxa\n")
        ntaxs = [10, 50, 100, 200, 500, 1000]
        hghts = [2000000]
        hghtlabs = ["medium"]
        rates = [0.000001]
        ratelabs = ["shallow"]
        
    elif do == "varyils":
        sys.stdout.write("Increasing ILS\n")
        ntaxs = [200]
        hghts = [10000000, 2000000, 500000]
        hghtlabs = ["low", "medium", "high"]
        rates = [0.0000001, 0.000001]
        ratelabs = ["deep", "shallow"]

    rows = []

    for ntax in ntaxs:
        for hghtlab, hght in zip(hghtlabs, hghts):
            for ratelab, rate in zip(ratelabs, rates):
                for ngen in ngens:
                    # Pick replicates
                    repls = [repl for repl in range(1, 51)]
                    if ntax == 10:
                        repls.remove(41)
                        if ngen == 1000:
                            repls.remove(12)  # additional because of concatenation
                    elif ntax == 50:
                        repls.remove(21)
                        repls.remove(27)  # additional because of concatenation
                        repls.remove(41)
                    elif ntax == 100:
                        repls.remove(8)
                        repls.remove(47)
                    elif ntax == 1000:
                        repls.remove(6)  # aster failure
                        repls.remove(8)  # aster failure
                        repls.remove(33) # aster failure
                        repls.remove(38) # aster failure
                    elif (hght == 500000) and (rate == 1e-6):
                        repls.remove(8)
                        repls.remove(15)
                        repls.remove(49)

                    # Pick methods
                    mthds = ["wastrid_s",
                             "wtreeqmc_wh_n2",
                             "aster_h",
                             "aster_h_t16"]

                    for repl in repls:
                        print("%d taxa, %d %g, %d genes : repl %d" % (ntax, hght, rate, ngen, repl))
                        for mthd in mthds:
                            for supp in supps:
                                xste_df = ste_df[(ste_df["NTAX"] == ntax) &
                                                 (ste_df["HGHT"] == hght) &
                                                 (ste_df["RATE"] == rate) &
                                                 (ste_df["REPL"] == repl) &
                                                 (ste_df["NGEN"] == ngen) &
                                                 (ste_df["MTHD"] == mthd) &
                                                 (ste_df["SUPP"] == supp)]

                                if xste_df.shape[0] != 1:
                                    print("  PROBLEM 1 - %s %s" % (mthd, supp))
                                    sefn = "NA"
                                    sefnr = "NA"
                                else:
                                    int1 = int(xste_df.I1.values[0])
                                    int2 = int(xste_df.I2.values[0])
                                    if int1 != int2:
                                        print("  POLY - %s %s!" % (mthd, supp))
                                    sefn = int(ste_df.FN.values[0])
                                    sefnr = float(sefn) / int1

                                xmrt_df = mrt_df[(mrt_df["NTAX"] == ntax) &
                                                 (mrt_df["HGHT"] == hght) &
                                                 (mrt_df["RATE"] == rate) &
                                                 (mrt_df["REPL"] == repl) &
                                                 (mrt_df["NGEN"] == ngen) &
                                                 (mrt_df["MTHD"] == mthd) &
                                                 (mrt_df["SUPP"] == supp)]

                                if xmrt_df.shape[0] != 1:
                                    print("  2 - PROBLEM - %s %s!" % (mthd, supp))
                                    node = "NA"
                                    secs = "NA"
                                else:
                                    node = xmrt_df.NODE.values[0]
                                    mrt = reformat_timing(xmrt_df.real.values[0])
                                    secs = mrt[0]

                                    if (mthd == mthds[0]):
                                        savenode = node
                                    else:
                                        if node != savenode:
                                            sys.exit("Methods run on different nodes!")

                                row = {}
                                row["NTAX"] = ntax
                                row["ILSL"] = hghtlab
                                row["SPEC"] = ratelab
                                row["REPL"] = repl
                                row["SUPP"] = supp
                                row["NGEN"] = ngen
                                row["MTHD"] = namemap[mthd]
                                row["NODE"] = node
                                row["SEFN"] = sefn
                                row["SEFNR"] = sefnr
                                row["SECS"] = secs
                                rows.append(row)

    ydf = pandas.DataFrame(rows, columns=cols)
    ydf.to_csv("data-" + do + "-error-and-runtime.csv",
               sep=',', na_rep="NA",header=True, index=False,
               quoting=csv.QUOTE_NONNUMERIC)
