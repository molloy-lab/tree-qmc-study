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


ste_df = pandas.read_csv("all_species_tree_error.csv", keep_default_na=False)
mrt_df = pandas.read_csv("all_runtime.csv", keep_default_na=False)
qrt_df = pandas.read_csv("all_quartets_runtime.csv", keep_default_na=False)
qsc_df = pandas.read_csv("all_quartet_score.csv", keep_default_na=False)

cols = ["NTAX", "STRHT", "SRATE", "REPL",
        "GTRE", "NGEN", "MTHD", "NODE",
        "SEFN", "SERF", "NQSC", "SECS", "FRAC"]

rows = []

ntaxs = [10, 50, 100]
strht = 2000000
srate = 0.000001
repls = range(1, 51)
ngens = [250, 1000]
gtres = ["estimatedgenetre", "truegenetrees"]

for ntax in ntaxs:
    for repl in repls:
        for ngen in ngens:
            for gtre in gtres:
                print("%d %d %d %s" % (ntax, repl, ngen, gtre))

                mthds = ["treeqmc_n0_v1.0.0", 
                         "treeqmc_n1_v1.0.0", 
                         "treeqmc_n2_v1.0.0",
                         "astral_3_v5.7.7", 
                         "fastral"]

                if gtre == "estimatedgenetre":
                    mthds.append("wqfm_v1.3")
                    mthds.append("wqmc_v3.0")

                    xqrt_df = qrt_df[(qrt_df["NTAX"] == ntax) &
                                     (qrt_df["STRHT"] == strht) &
                                     (qrt_df["SRATE"] == srate) &
                                     (qrt_df["REPL"] == repl) & 
                                     (qrt_df["GTRE"] == gtre) &
                                     (qrt_df["NGEN"] == ngen)]

                    if xqrt_df.shape[0] != 1:
                        print(xqrt_df)
                        sys.exit("1 ERROR!")

                    qrt = reformat_timing(xqrt_df.real.values[0])

                for mthd in mthds:
                    xste_df = ste_df[(ste_df["NTAX"] == ntax) &
                                     (ste_df["STRHT"] == strht) &
                                     (ste_df["SRATE"] == srate) &
                                     (ste_df["REPL"] == repl) & 
                                     (ste_df["GTRE"] == gtre) &
                                     (ste_df["NGEN"] == ngen) &
                                     (ste_df["MTHD"] == mthd)]

                    if xste_df.shape[0] != 1:
                        sys.exit("2 ERROR!")

                    serf = float(xste_df.RF.values[0])
                    sefn = int(xste_df.FN.values[0])

                    xmrt_df = mrt_df[(mrt_df["NTAX"] == ntax) &
                                     (mrt_df["STRHT"] == strht) &
                                     (mrt_df["SRATE"] == srate) &
                                     (mrt_df["REPL"] == repl) &
                                     (mrt_df["GTRE"] == gtre) &
                                     (mrt_df["NGEN"] == ngen) &
                                     (mrt_df["MTHD"] == mthd)]

                    if xmrt_df.shape[0] != 1:
                        sys.exit("3 ERROR!")

                    node = xmrt_df.NODE.values[0]
                    mrt = reformat_timing(xmrt_df.real.values[0])
                    secs = mrt[0]

                    if (mthd == "treeqmc_n0_v1.0.0") and (gtre == "estimatedgenetre"):
                        savenode = node
                    else:
                        if node != savenode:
                            sys.exit("Methods run on different nodes!")

                    xqsc_df = qsc_df[(qsc_df["NTAX"] == ntax) &
                                     (qsc_df["STRHT"] == strht) &
                                     (qsc_df["SRATE"] == srate) &
                                     (qsc_df["REPL"] == repl) & 
                                     (qsc_df["GTRE"] == gtre) &
                                     (qsc_df["NGEN"] == ngen) &
                                     (qsc_df["MTHD"] == mthd)]

                    if xqsc_df.shape[0] != 1:
                        sys.exit("4 ERROR!")

                    nqsc = float(xqsc_df.NQS.values[0])

                    row = {}
                    row["NTAX"] = ntax
                    row["STRHT"] = strht
                    row["SRATE"] = srate
                    row["REPL"] = repl
                    row["GTRE"] = gtre
                    row["NGEN"] = ngen
                    row["MTHD"] = mthd
                    row["NODE"] = node
                    row["SEFN"] = sefn
                    row["SERF"] = serf
                    row["NQSC"] = nqsc
                    row["SECS"] = secs
                    row["FRAC"] = 0
                    if (mthd == "wqfm_v1.3") or (mthd == "wqmc_v3.0"):
                        row["SECS"] += qrt[0]
                        row["FRAC"] = qrt[0] / row["SECS"]  # Fraction of time spent on quartets
                    rows.append(row)

                # Now append row for true species tree
                xqsc_df = qsc_df[(qsc_df["NTAX"] == ntax) &
                                 (qsc_df["STRHT"] == strht) &
                                 (qsc_df["SRATE"] == srate) &
                                 (qsc_df["REPL"] == repl) & 
                                 (qsc_df["GTRE"] == gtre) &
                                 (qsc_df["NGEN"] == ngen) &
                                 (qsc_df["MTHD"] == "TRUE")]

                if xqsc_df.shape[0] != 1:
                    sys.exit("5 ERROR!")

                nqsc = float(xqsc_df.NQS.values[0])

                row = {}
                row["NTAX"] = ntax
                row["STRHT"] = strht
                row["SRATE"] = srate
                row["REPL"] = repl
                row["GTRE"] = gtre
                row["NGEN"] = ngen
                row["MTHD"] = "TRUE"
                row["NODE"] = node
                row["SEFN"] = 0
                row["SERF"] = 0
                row["NQSC"] = nqsc
                row["SECS"] = 0
                row["FRAC"] = 0
                rows.append(row)


df = pandas.DataFrame(rows, columns=cols)
df.to_csv("data-lt200tax-error-and-timings.csv",
          sep=',', na_rep="NA",header=True, index=False)
