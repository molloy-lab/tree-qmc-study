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
wqf_df = pandas.read_csv("all_wqfit.csv", keep_default_na=False)
mrt_df = pandas.read_csv("all_runtime.csv", keep_default_na=False)
qrt_df = pandas.read_csv("all_quartets_runtime.csv", keep_default_na=False)
qsc_df = pandas.read_csv("all_quartet_score.csv", keep_default_na=False)

cols = ["SCAL", "NGEN", "NBPS", "REPL",
        "MTHD", "NODE", "SEFN", "SERF", 
        "WQFT", "NQSC", "SECS", "FRAC"]

rows = []

repls = [str("R%d" % x) for x in range(1, 21)]

for do in ["ils", "ngen"]:
    if do == "ils":
        scals =["0.5X", "1X", "2X"]
        ngens = [1000]
        nbpss = ["true", "500"]
    elif do == "ngen":
        scals = ["1X"]
        ngens = [50, 100, 200, 500]
        nbpss = ["true", "500"]

    for scal in scals:
        for nbps in nbpss:
            for ngen in ngens:
                for repl in repls:
                    print("%s %s %s %s" % (scal, ngen, nbps, repl))

                    mthds = ["treeqmc_n0_v1.0.0", 
                             "treeqmc_n1_v1.0.0", 
                             "treeqmc_n2_v1.0.0",
                             "astral_3_v5.7.7", 
                             "fastral",
                             "wqfm_v1.3",
                             "wqmc_v3.0"]

                    # Process quartet runtime
                    xqrt_df = qrt_df[(qrt_df["SCAL"] == scal) &
                                     (qrt_df["NGEN"] == ngen) &
                                     (qrt_df["NBPS"] == nbps) &
                                     (qrt_df["REPL"] == repl)]

                    if xqrt_df.shape[0] != 1:
                        sys.exit("1 ERROR - %s qrt!" % mthd)

                    qrt = reformat_timing(xqrt_df.real.values[0])

                    for mthd in mthds:
                        # Process species tree error
                        xste_df = ste_df[(ste_df["SCAL"] == scal) &
                                         (ste_df["NGEN"] == ngen) &
                                         (ste_df["NBPS"] == nbps) &
                                         (ste_df["REPL"] == repl) &
                                         (ste_df["MTHD"] == mthd)]

                        if xste_df.shape[0] != 1:
                            sys.exit("2 ERROR - %s ste!" % mthd)

                        serf = float(xste_df.RF.values[0])
                        sefn = int(xste_df.FN.values[0])

                        # Process wQfit
                        xwqf_df = wqf_df[(wqf_df["SCAL"] == scal) &
                                         (wqf_df["NGEN"] == ngen) &
                                         (wqf_df["NBPS"] == nbps) &
                                         (wqf_df["REPL"] == repl) &
                                         (wqf_df["MTHD"] == mthd)]

                        if xwqf_df.shape[0] != 1:
                            sys.exit("3 ERROR - %s ste!" % mthd)

                        wqft = float(xwqf_df.WQFIT.values[0])

                        # Process runtime
                        xmrt_df = mrt_df[(mrt_df["SCAL"] == scal) &
                                         (mrt_df["NGEN"] == ngen) &
                                         (mrt_df["NBPS"] == nbps) &
                                         (mrt_df["REPL"] == repl) &
                                         (mrt_df["MTHD"] == mthd)]

                        if xmrt_df.shape[0] != 1:
                            sys.exit("4 ERROR - %s rt!" % mthd)

                        node = xmrt_df.NODE.values[0]
                        mrt = reformat_timing(xmrt_df.real.values[0])
                        secs = mrt[0]

                        if (mthd == "treeqmc_n0_v1.0.0"):
                            savenode = node
                        else:
                            if node != savenode:
                                sys.stdout.write("Methods run on different nodes!")

                        # Process quartet score
                        xqsc_df = qsc_df[(qsc_df["SCAL"] == scal) &
                                         (qsc_df["NGEN"] == ngen) &
                                         (qsc_df["NBPS"] == nbps) &
                                         (qsc_df["REPL"] == repl) &
                                         (qsc_df["MTHD"] == mthd)]

                        if xqsc_df.shape[0] != 1:
                            sys.exit("5 ERROR - %s qsc!" % mthd)

                        nqsc = float(xqsc_df.NQS.values[0])

                        row = {}
                        row["SCAL"] = scal
                        row["NGEN"] = ngen
                        row["NBPS"] = nbps
                        row["REPL"] = repl
                        row["MTHD"] = mthd
                        row["NODE"] = node
                        row["SEFN"] = sefn
                        row["WQFT"] = wqft
                        row["SERF"] = serf
                        row["NQSC"] = nqsc
                        row["SECS"] = secs
                        row["FRAC"] = 0
                        if (mthd == "wqfm_v1.3") or (mthd == "wqmc_v3.0"):
                            row["SECS"] += qrt[0]
                            row["FRAC"] = qrt[0] / row["SECS"]  # Fraction of time spent on quartets
                        rows.append(row)

                    # Now append row for true species tree
                    xqsc_df = qsc_df[(qsc_df["SCAL"] == scal) &
                                     (qsc_df["NGEN"] == ngen) &
                                     (qsc_df["NBPS"] == nbps) &
                                     (qsc_df["REPL"] == repl) & 
                                     (qsc_df["MTHD"] == "TRUE")]

                    if xqsc_df.shape[0] != 1:
                        sys.exit("6 ERROR - TRUE!")

                    nqsc = float(xqsc_df.NQS.values[0])

                    row = {}
                    row["SCAL"] = scal
                    row["NGEN"] = ngen
                    row["NBPS"] = nbps
                    row["REPL"] = repl
                    row["MTHD"] = "TRUE"
                    row["NODE"] = node
                    row["SEFN"] = 0
                    row["SERF"] = 0
                    row["WQFT"] = 1
                    row["NQSC"] = nqsc
                    row["SECS"] = 0
                    row["FRAC"] = 0
                    rows.append(row)

df = pandas.DataFrame(rows, columns=cols)
df.to_csv("data-all-error-and-timings.csv",
          sep=',', na_rep="NA",header=True, index=False)

