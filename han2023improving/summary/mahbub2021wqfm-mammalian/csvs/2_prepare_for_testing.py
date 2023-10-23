import pandas
import math
import numpy
import sys


df = pandas.read_csv("data-all-error-and-timings.csv")
df.fillna("NA", inplace=True)

cols = ["SCAL", "NGEN", "NBPS", "REPL",
        "NODE",
        "wQFMxSEFN", "wQFMxSERF", "wQFMxNQSC", "wQFMxSECS",
        "TQMCn2xSEFN",  "TQMCn2xSERF",  "TQMCn2xNQSC",  "TQMCn2xSECS",
        "FASTRALxSEFN", "FASTRALxSERF", "FASTRALxNQSC", "FASTRALxSECS",
        "ASTRAL3xSEFN", "ASTRAL3xSERF", "ASTRAL3xNQSC", "ASTRAL3xSECS"]
rows = []

scals = []
ngens = []
nbpss = []
mthds = ["wqfm_v1.3",
         "treeqmc_n2_v1.0.0",
         "astral_3_v5.7.7",
         "fastral"]

for scal in ["noscale", "scale2d", "scale2u"]:
    if scal == "noscale":
        ngens = ["25g", "50g", "100g", "200g", "400g", "800g"]
    else:
        ngens = ["200g"]
    for ngen in ngens:
        if (scal == "noscale") and (ngen == "200g"):
            nbpss = ["250b", "500b", "1000b", "1500b", "true"]
        else:
             nbpss = ["500b"]

        for nbps in nbpss:
            repls = df[(df["SCAL"] == scal) &
                       (df["NGEN"] == ngen) &
                       (df["NBPS"] == nbps) &
                       (df["MTHD"] == "treeqmc_n2_v1.0.0")].REPL.values
            for repl in repls:
                print("%s %s %s" % (scal, ngen, nbps))
                xdf = df[(df["SCAL"] == scal) &
                         (df["NGEN"] == ngen) &
                         (df["NBPS"] == nbps) &
                         (df["REPL"] == repl)]

                wqfm = xdf[(xdf["MTHD"] == "wqfm_v1.3")]
                astral = xdf[(xdf["MTHD"] == "astral_3_v5.7.7")]
                fastral = xdf[(xdf["MTHD"] == "fastral")]
                tqmc2 = xdf[(xdf["MTHD"] == "treeqmc_n2_v1.0.0")]

                node1 = tqmc2.NODE.values[0]
                for mthd in mthds:
                    node2 = xdf[(xdf["MTHD"] == mthd)].NODE.values[0]
                    if (node1 != node2) and (node2 != "NA"):
                        print(xdf)
                        sys.exit("Node error!")

                row = {}
                row["SCAL"] = scal
                row["NGEN"] = ngen
                row["NBPS"] = nbps
                row["REPL"] = repl
                row["NODE"] = node1

                row["wQFMxSEFN"] = wqfm.SEFN.values[0]
                row["wQFMxSERF"] = wqfm.SERF.values[0]
                row["wQFMxNQSC"] = wqfm.NQSC.values[0]
                row["wQFMxSECS"] = wqfm.SECS.values[0]

                row["TQMCn2xSEFN"] = tqmc2.SEFN.values[0]
                row["TQMCn2xSERF"] = tqmc2.SERF.values[0]
                row["TQMCn2xNQSC"] = tqmc2.NQSC.values[0]
                row["TQMCn2xSECS"] = tqmc2.SECS.values[0]

                row["FASTRALxSEFN"] = fastral.SEFN.values[0]
                row["FASTRALxSERF"] = fastral.SERF.values[0]
                row["FASTRALxNQSC"] = fastral.NQSC.values[0]
                row["FASTRALxSECS"] = fastral.SECS.values[0]

                row["ASTRAL3xSEFN"] = astral.SEFN.values[0]
                row["ASTRAL3xSERF"] = astral.SERF.values[0]
                row["ASTRAL3xNQSC"] = astral.NQSC.values[0]
                row["ASTRAL3xSECS"] = astral.SECS.values[0]

                rows.append(row)

df = pandas.DataFrame(rows, columns=cols)

df.to_csv("data-for-testing.csv",
          sep=',', na_rep="NA",header=True, index=False)

