import pandas
import math
import numpy
import sys

#sys.exit("DONE RUNNING")

df = pandas.read_csv("data-all-error.csv")
df.fillna("NA", inplace=True)

cols = ["NBPS", "NGEN", "SUPP", "REPL", 
        "WASTRIDxSEFN", "WASTRIDxSERF",
        "ASTERHxSEFN", "ASTERHxSERF",
        "TQMCwhn2xSEFN", "TQMCwhn2xSERF",
        "TQMCwhn1xSEFN", "TQMCwhn1xSERF",
        "TQMCwhn0xSEFN", "TQMCwhn0xSERF",
        "TQMCwsn2xSEFN", "TQMCwsn2xSERF",
        "TQMCwnn2xSEFN", "TQMCwnn2xSERF"]

rows = []

mthds = ["wastrid_s",
         "aster_h_t16",
         "wtreeqmc_wh_n2",
         "wtreeqmc_wh_n1",
         "wtreeqmc_wh_n0",
         "wtreeqmc_ws_n2",
         "wtreeqmc_wn_n2"]

nbpss = [200, 400, 800, 1600]
ngens = [50, 200, 500, 1000]
supps = ["bs", "abayes"]

for nbps in nbpss:
    for ngen in ngens:
        for supp in supps:
            repls = df[(df["NBPS"] == nbps) &
                       (df["NGEN"] == ngen) & 
                       (df["SUPP"] == supps[0]) & 
                       (df["MTHD"] == "wtreeqmc_wh_n2")].REPL.values

            for repl in repls:
                print("%d %d %s %d" % (nbps, ngen, supp, repl))

                xdf = df[(df["NBPS"] == nbps) &
                         (df["NGEN"] == ngen) & 
                         (df["REPL"] == repl)]

                wastrid = xdf[(xdf["MTHD"] == "wastrid_s") & (xdf["SUPP"] == supp)]
                asterh = xdf[(xdf["MTHD"] == "aster_h_t16") & (xdf["SUPP"] == supp)]
                tqmc_whn2 = xdf[(xdf["MTHD"] == "wtreeqmc_wh_n2") & (xdf["SUPP"] == supp)]
                tqmc_whn1 = xdf[(xdf["MTHD"] == "wtreeqmc_wh_n1") & (xdf["SUPP"] == supp)]
                tqmc_whn0 = xdf[(xdf["MTHD"] == "wtreeqmc_wh_n0") & (xdf["SUPP"] == supp)]
                tqmc_wsn2 = xdf[(xdf["MTHD"] == "wtreeqmc_ws_n2") & (xdf["SUPP"] == supp)]
                tqmc_wnn2 = xdf[(xdf["MTHD"] == "wtreeqmc_wn_n2") & (xdf["SUPP"] == "none_refinepoly")]  ## LIKE ORIGINAL

                row = {}
                row["NBPS"] = nbps
                row["NGEN"] = ngen
                row["SUPP"] = supp
                row["REPL"] = repl

                row["WASTRIDxSEFN"] = wastrid.SEFN.values[0]
                row["WASTRIDxSERF"] = wastrid.SERF.values[0]

                row["ASTERHxSEFN"] = asterh.SEFN.values[0]
                row["ASTERHxSERF"] = asterh.SERF.values[0]

                row["TQMCwhn2xSEFN"] = tqmc_whn2.SEFN.values[0]
                row["TQMCwhn2xSERF"] = tqmc_whn2.SERF.values[0]

                row["TQMCwhn1xSEFN"] = tqmc_whn1.SEFN.values[0]
                row["TQMCwhn1xSERF"] = tqmc_whn1.SERF.values[0]

                row["TQMCwhn0xSEFN"] = tqmc_whn0.SEFN.values[0]
                row["TQMCwhn0xSERF"] = tqmc_whn0.SERF.values[0]

                row["TQMCwsn2xSEFN"] = tqmc_wsn2.SEFN.values[0]
                row["TQMCwsn2xSERF"] = tqmc_wsn2.SERF.values[0]

                row["TQMCwnn2xSEFN"] = tqmc_wnn2.SEFN.values[0]
                row["TQMCwnn2xSERF"] = tqmc_wnn2.SERF.values[0]

                rows.append(row)

df = pandas.DataFrame(rows, columns=cols)

df.to_csv("data-for-testing.csv",
          sep=',', na_rep="NA",header=True, index=False)
