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
+ 3 remove 1 replicate (i.e. 41) from 10-taxon
+ 6 remove 2 replicates (i.e. 21, 41) from 50-taxon (note 27 was close with 590ish genes)
+ 6 remove 2 replicates (i.e. 8 and 47) from 100-taxon
+ 9 remove 3 replicates (i.e. 8, 15, 49) from 200-taxon/500K/1e-06

Concatenation (CA-ML) tree was also missing on 
+ 3 remove 1 replicate (i.e. 27) for 50-taxon, related to above
+ 1 remove 1 replicate (i.e. 12) for 10-taxon with 1000 genes
+ Also CA-ML wasn't run on any 1000-taxon with 1000 genes

TREE-QMC returned polytomies on 
+ replicate 41 for 10-taxon with 50 genes (but already removed)
+ replicate 8 for 200-taxon/500K/1e-06 with 200 genes (but already removed)

Failures due to wall clock time for single-threaded ASTER-h:
+ 2 remove 2 replicates (8 and 38) from 1000-taxon with 1000 genes
+ 1 remove 1 replicate (21) from 200-taxon/500K/1e-07 with 1000 genes
"""

sys.exit("DONE RUNNING")

namemap = {}
namemap["aster_h_t16"] = "ASTER-wh (16 threads)"
namemap["aster_h"] = "ASTER-wh (1 thread)"
namemap["wastrid_s"] = "ASTRID-ws"
namemap["asteroid"] = "ASTEROID"
namemap["wtreeqmc_wf_n2"] = "TQMC-n2"
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

mrt_df = pandas.read_csv("all_runtime.csv.gz", keep_default_na=False, compression="gzip")

cols = ["NTAX", "ILSL", "SPEC", "REPL", "SUPP", "NGEN", "MTHD", "NODE", "SECS"]

ngens = [50, 200, 1000]

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
                        repls.remove(41)      # +3
                        if ngen == 1000:
                            repls.remove(12)  # +1 additional because of concatenation
                    elif ntax == 50:
                        repls.remove(21)  # +3
                        repls.remove(27)  # +3 additional because of concatenation
                        repls.remove(41)  # +3
                    elif ntax == 100:
                        repls.remove(8)   # +3
                        repls.remove(47)  # +3
                    elif (ntax == 1000) and (ngen == 1000):
                            repls.remove(8)  # +1 additional because of aster single threaded failure
                            repls.remove(38) # +1 additional because of aster single threaded failure
                    elif (hght == 500000):
                        if (rate == 1e-6):
                            repls.remove(8)  # +3
                            repls.remove(15) # +3
                            repls.remove(49) # +3
                        else:
                            if ngen == 1000:
                                repls.remove(21) # +1 additional because of aster single threaded failure

                    # Pick methods
                    mthds = ["wastrid_s",
                             "asteroid",
                             "wtreeqmc_wh_n2",
                             "aster_h_t16",
                             "aster_h",
                             "wtreeqmc_wf_n2"] # ran separately

                    for repl in repls:
                        print("%d taxa, %d %g, %d genes : repl %d" % (ntax, hght, rate, ngen, repl))
                        total = 0
                        for mthd in mthds:
                            if mthd == "wtreeqmc_wf_n2":
                                supps = ["sh", "abayes"]
                            else:
                                supps = ["abayes"]

                            for supp in supps:
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

                                    if (mthd != "aster_h") and (mthd != "wtreeqmc_wf_n2"):
                                        total += secs

                                    #if (mthd == mthds[0]):
                                    #    savenode = node
                                    #else:
                                    #    if node != savenode:
                                    #        sys.exit("Methods run on different nodes!")

                                

                                row = {}
                                row["NTAX"] = ntax
                                row["ILSL"] = hghtlab
                                row["SPEC"] = ratelab
                                row["REPL"] = repl
                                row["SUPP"] = supp
                                row["NGEN"] = ngen
                                row["MTHD"] = namemap[mthd]
                                row["NODE"] = node
                                row["SECS"] = secs
                                rows.append(row)

                        if (total / (60 * 60)) > 4:
                            print(total)
                            sys.exit("WARNING: Checking ASTER got at least 20 hours")

    ydf = pandas.DataFrame(rows, columns=cols)
    ydf.to_csv("data-" + do + "-runtime.csv",
               sep=',', na_rep="NA",header=True, index=False,
               quoting=csv.QUOTE_NONNUMERIC)
