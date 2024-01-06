import pandas
import math
import numpy
import sys

sys.exit("DONE RUNNING")

cols = ["NTAX", "ILSL", "SPEC", "REPL", "SUPP", "NGEN",
        "CAMLxSEFN", "CAMLxSERF",
        "WASTRIDxSEFN", "WASTRIDxSERF",
        "ASTERHxSEFN", "ASTERHxSERF",
        "TQMCwhn2xSEFN", "TQMCwhn2xSERF",
        "TQMCwsn2xSEFN", "TQMCwsn2xSERF",
        "TQMCwnn2xSEFN", "TQMCwnn2xSERF"]

ngens = [50, 200, 1000]
supps = ["sh", "abayes"]

for do in ["varyntax", "varyils"]: 
    if do == "varyntax":
        sys.stdout.write("Increasing number of taxa\n")
        ntaxs = [10, 50, 100, 200, 500, 1000]
        hghts = ["medium"]
        rates = ["shallow"]
        
    elif do == "varyils":
        sys.stdout.write("Increasing ILS\n")
        ntaxs = [200]
        hghts = ["low", "medium", "high"]
        rates = ["deep", "shallow"]

    df = pandas.read_csv("data-" + do + "-error.csv", keep_default_na=False)
    df.fillna("NA", inplace=True)
    rows = []

    for ntax in ntaxs:
        for hght in hghts:
            for rate in rates:
                for ngen in ngens:
                    xdf = df[(df["NTAX"] == ntax) &
                             (df["ILSL"] == hght) &
                             (df["SPEC"] == rate) &
                             (df["NGEN"] == ngen)]

                    repls = xdf[(xdf["MTHD"] == "CA-ML")].REPL.values

                    for repl in repls:
                        ydf = xdf[(xdf["REPL"] == repl)]

                        caml = ydf[(ydf["MTHD"] == "CA-ML")]
                        tqmc_wnn2 = ydf[(ydf["MTHD"] == "TQMC-n2")]

                        for supp in supps:
                            print("%d %s %s %d %d %s" % (ntax, hght, rate, repl, ngen, supp))
                            zdf = ydf[(ydf["SUPP"] == supp)]

                            wastrid = zdf[(zdf["MTHD"] == "ASTRID-ws")]
                            asterh = zdf[(zdf["MTHD"] == "ASTER-wh")]
                            tqmc_whn2 = zdf[(zdf["MTHD"] == "TQMC-wh_n2")]
                            tqmc_wsn2 = zdf[(zdf["MTHD"] == "TQMC-ws_n2")]

                            row = {}
                            row["NTAX"] = ntax
                            row["ILSL"] = hght
                            row["SPEC"] = rate
                            row["REPL"] = repl
                            row["SUPP"] = supp
                            row["NGEN"] = ngen

                            row["CAMLxSEFN"] = caml.SEFN.values[0]
                            row["CAMLxSERF"] = caml.SEFNR.values[0]

                            row["WASTRIDxSEFN"] = wastrid.SEFN.values[0]
                            row["WASTRIDxSERF"] = wastrid.SEFNR.values[0]

                            row["ASTERHxSEFN"] = asterh.SEFN.values[0]
                            row["ASTERHxSERF"] = asterh.SEFNR.values[0]

                            row["TQMCwhn2xSEFN"] = tqmc_whn2.SEFN.values[0]
                            row["TQMCwhn2xSERF"] = tqmc_whn2.SEFNR.values[0]

                            row["TQMCwsn2xSEFN"] = tqmc_wsn2.SEFN.values[0]
                            row["TQMCwsn2xSERF"] = tqmc_wsn2.SEFNR.values[0]

                            row["TQMCwnn2xSEFN"] = tqmc_wnn2.SEFN.values[0]
                            row["TQMCwnn2xSERF"] = tqmc_wnn2.SEFNR.values[0]

                            rows.append(row)

    df = pandas.DataFrame(rows, columns=cols)

    df.to_csv("data-" + do + "-for-testing.csv",
            sep=',', na_rep="NA",header=True, index=False)
