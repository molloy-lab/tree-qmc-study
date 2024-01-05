import pandas
import math
import numpy
import sys

sys.exit("DONE RUNNING")

cols = ["NTAX", "HGHT", "RATE", "REPL", "SUPP", "NGEN",
        "CAMLxSEFN", "CAMLxSERF",
        "WASTRIDxSEFN", "WASTRIDxSERF",
        "ASTERHxSEFN", "ASTERHxSERF",
        "TQMCwhn2xSEFN", "TQMCwhn2xSERF",
        "TQMCwhn1xSEFN", "TQMCwhn1xSERF",
        "TQMCwhn0xSEFN", "TQMCwhn0xSERF",
        "TQMCwsn2xSEFN", "TQMCwsn2xSERF",
        "TQMCwnn2xSEFN", "TQMCwnn2xSERF"]
rows = []

ngens = [50, 200, 1000]
supps = ["sh", "abayes"]
mthds = ["caml",
         "asteroid",
         "wastrid_s",
         "aster_h_t16",
         "wtreeqmc_wh_n2",
         "wtreeqmc_wh_n1",
         "wtreeqmc_wh_n0",
         "wtreeqmc_ws_n2",
         "wtreeqmc_wn_n2"]

for do in ["ntax", "ils"]: 
    if do == "ntax":
        sys.stdout.write("Increasing number of taxa\n")
        df = pandas.read_csv("data-varyntax-error.csv", keep_default_na=False)
        ntaxs = [10, 50, 100, 500, 1000]
        hghts = [2000000]
        rates = [0.000001]
    elif do == "ils":
        sys.stdout.write("Increasing ILS\n")
        df = pandas.read_csv("data-varyils-error.csv", keep_default_na=False)
        ntaxs = [200]
        hghts = [10000000, 2000000, 500000]
        rates = [0.0000001, 0.000001]

    df.fillna("NA", inplace=True)

    for ntax in ntaxs:
        for hght in hghts:
            for rate in rates:
                for ngen in ngens:
                    xdf = df[(df["NTAX"] == ntax) &
                             (df["HGHT"] == hght) &
                             (df["RATE"] == rate) &
                             (df["NGEN"] == ngen)]

                    repls = xdf[(xdf["MTHD"] == "caml")].REPL.values

                    for repl in repls:
                        ydf = xdf[(xdf["REPL"] == repl)]

                        caml = ydf[(ydf["MTHD"] == "caml")]
                        tqmc_wnn2 = ydf[(ydf["MTHD"] == "wtreeqmc_wn_n2")]

                        for supp in supps:
                            print("%d %d %g %d %d %s" % (ntax, hght, rate, repl, ngen, supp))
                            zdf = ydf[(ydf["SUPP"] == supp)]

                            wastrid = zdf[(zdf["MTHD"] == "wastrid_s")]
                            asterh = zdf[(zdf["MTHD"] == "aster_h_t16")]
                            tqmc_whn2 = zdf[(zdf["MTHD"] == "wtreeqmc_wh_n2")]
                            tqmc_whn1 = zdf[(zdf["MTHD"] == "wtreeqmc_wh_n1")]
                            tqmc_whn0 = zdf[(zdf["MTHD"] == "wtreeqmc_wh_n0")]
                            tqmc_wsn2 = zdf[(zdf["MTHD"] == "wtreeqmc_ws_n2")]

                            row = {}
                            row["NTAX"] = ntax
                            row["HGHT"] = hght
                            row["RATE"] = rate
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

                            row["TQMCwhn1xSEFN"] = tqmc_whn1.SEFN.values[0]
                            row["TQMCwhn1xSERF"] = tqmc_whn1.SEFNR.values[0]

                            row["TQMCwhn0xSEFN"] = tqmc_whn0.SEFN.values[0]
                            row["TQMCwhn0xSERF"] = tqmc_whn0.SEFNR.values[0]

                            row["TQMCwsn2xSEFN"] = tqmc_wsn2.SEFN.values[0]
                            row["TQMCwsn2xSERF"] = tqmc_wsn2.SEFNR.values[0]

                            row["TQMCwnn2xSEFN"] = tqmc_wnn2.SEFN.values[0]
                            row["TQMCwnn2xSERF"] = tqmc_wnn2.SEFNR.values[0]

                            rows.append(row)

df = pandas.DataFrame(rows, columns=cols)

df.to_csv("data-for-testing.csv",
          sep=',', na_rep="NA",header=True, index=False)
