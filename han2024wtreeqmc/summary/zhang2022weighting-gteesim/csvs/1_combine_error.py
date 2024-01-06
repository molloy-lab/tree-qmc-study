import pandas
import numpy
import sys

sys.exit("DONE RUNNING")

ste_df = pandas.read_csv("all_species_tree_error.csv", keep_default_na=False)

cols = ["NBPS", "NGEN", "REPL", "SUPP",
        "MTHD", "SEFN", "SEFP", "SERF"]

rows = []

mthds = ["aster_h_t16",
         "wastrid_s",
         "wtreeqmc_wh_n2",
         "wtreeqmc_wn_n2",
         "wtreeqmc_ws_n2",
         "wtreeqmc_wh_n1",
         "wtreeqmc_wh_n0"]

namemap = {}
namemap["aster_h_t16"] = "ASTER-wh"
namemap["wastrid_s"] = "ASTRID-ws"
namemap["wtreeqmc_wn_n2"] = "TQMC-wn_n2"
namemap["wtreeqmc_ws_n2"] = "TQMC-ws_n2"
namemap["wtreeqmc_wh_n2"] = "TQMC-wh_n2"
namemap["wtreeqmc_wh_n1"] = "TQMC-wh_n1"
namemap["wtreeqmc_wh_n0"] = "TQMC-wh_n0"

nbpss = [200, 400, 800, 1600]
ngens = [50, 200, 500, 1000]
repls = range(1, 51)
supps = ["bs", "abayes"]

for nbps in nbpss:
    for ngen in ngens:
        for repl in repls:
            print("%d %d %d" % (nbps, ngen, repl))
            for mthd in mthds:
                for supp in supps:
                    xste_df = ste_df[(ste_df["NBPS"] == nbps) &
                                     (ste_df["NGEN"] == ngen) &
                                     (ste_df["REPL"] == repl) &
                                     (ste_df["SUPP"] == supp) &
                                     (ste_df["MTHD"] == mthd)]

                    if xste_df.shape[0] != 1:
                        sys.exit("  1 ERROR - %s!\n" % mthd)

                    if (xste_df.FN.values[0] == '') or (xste_df.FN.values[0] == "NA"):
                        sefn = "NA"
                        sefp = "NA"
                        serf = "NA"
                    else:
                        int1 = int(xste_df.I1.values[0])
                        int2 = int(xste_df.I2.values[0])
                        if int1 != int2:
                            sys.exit("  2 ERROR - %s!\n" % mthd)
                        sefn = int(xste_df.FN.values[0])
                        sefp = int(xste_df.FP.values[0])
                        serf = float(sefn) / int1

                    row = {}
                    row["NBPS"] = nbps
                    row["NGEN"] = ngen
                    row["REPL"] = repl
                    if mthd == "wtreeqmc_wn_n2":
                        if supp == "bs":
                            row["SUPP"] = "none_keeppoly"
                        else:
                            row["SUPP"] = "none_refinepoly"
                    else:
                        row["SUPP"] = supp
                    row["MTHD"] = namemap[mthd]
                    row["SEFN"] = sefn
                    row["SEFP"] = sefp
                    row["SERF"] = serf
                    rows.append(row)

df = pandas.DataFrame(rows, columns=cols)
df.to_csv("data-all-error.csv",
          sep=',', na_rep="NA",header=True, index=False)
