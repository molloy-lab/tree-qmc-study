import pandas
import numpy
import sys


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
qsc_df = pandas.read_csv("all_quartet_score.csv", keep_default_na=False)

cols = ["NTAX", "STRHT", "SRATE", "REPL",
        "GTRE", "NGEN", "MTHD", "NODE",
        "SEFN", "SERF", "NQSC", "SECS", "FRAC"]

rows = []

ntaxs = [500, 1000]
strht = 2000000
srate = 0.000001
repls = range(1, 51)
gtres = ["estimatedgenetre", "truegenetrees"]
ngens = [250, 1000]
mthds = ["treeqmc_n0_v1.0.0",
         "treeqmc_n1_v1.0.0",
         "treeqmc_n2_v1.0.0",
         "astral_3_v5.7.7",
         "fastral"]


for ntax in ntaxs:
    for repl in repls:
        for ngen in ngens:
            for gtre in gtres:
                print("%d %d %s %d" % (ntax, repl, gtre, ngen))
                for mthd in mthds:
                    xste_df = ste_df[(ste_df["NTAX"] == ntax) &
                                     (ste_df["STRHT"] == strht) &
                                     (ste_df["SRATE"] == srate) &
                                     (ste_df["REPL"] == repl) & 
                                     (ste_df["GTRE"] == gtre) &
                                     (ste_df["NGEN"] == ngen) &
                                     (ste_df["MTHD"] == mthd)]

                    if xste_df.shape[0] != 1:
                        #sys.exit("1 ERROR!")
                        print(" - 1 bad %s" % mthd)
                        serf = "NA"
                        sefn = "NA"
                    else:
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
                        #sys.exit("2 ERROR!")
                        print(" - 2 bad %s" % mthd)
                        node = "NA"
                        secs = "NA"
                    else:
                        node = xmrt_df.NODE.values[0]
                        mrt = reformat_timing(xmrt_df.real.values[0])
                        secs = mrt[0]

                    if (mthd == "treeqmc_n0_v1.0.0") and (gtre == "estimatedgenetre"):
                        savenode = node
                    else:
                        if node != savenode:
                            sys.stdout.write(" - Method %s run on different node!\n" % mthd)

                    xqsc_df = qsc_df[(qsc_df["NTAX"] == ntax) &
                                     (qsc_df["STRHT"] == strht) &
                                     (qsc_df["SRATE"] == srate) &
                                     (qsc_df["REPL"] == repl) & 
                                     (qsc_df["GTRE"] == gtre) &
                                     (qsc_df["NGEN"] == ngen) &
                                     (qsc_df["MTHD"] == mthd)]

                    if xqsc_df.shape[0] != 1:
                        #sys.exit("3 ERROR!")
                        print(" - 3 bad %s" % mthd)
                        nqsc = "NA"
                    else:
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
                    row["FRAC"] = 0.0
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
df.to_csv("data-gt200tax-error-and-timings.csv",
          sep=',', na_rep="NA", header=True, index=False)

