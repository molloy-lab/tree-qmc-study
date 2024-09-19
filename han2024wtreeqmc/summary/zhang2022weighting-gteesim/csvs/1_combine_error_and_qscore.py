import pandas
import numpy
import sys

sys.exit("DONE RUNNING")

ste_df = pandas.read_csv("all_species_tree_error.csv.gz", keep_default_na=False, compression="gzip")
qsc_df = pandas.read_csv("all_quartet_score.csv.gz", keep_default_na=False, compression="gzip")

cols = ["NBPS", "NGEN", "REPL", "SUPP",
        "MTHD", "SEFN", "SEFP", "SERF",
        "QSWN", "QSWH", "TRUE_QSWN", "TRUE_QSWH"]
rows = []

mthds = ["aster_h",
         "wastrid_s",
         "wtreeqmc_wf_n2",
         "wtreeqmc_wn_n2",
         "wtreeqmc_ws_n2",
         "wtreeqmc_wh_n2",
         "wtreeqmc_wh_n1",
         "wtreeqmc_wh_n0"]

namemap = {}
namemap["aster_h"] = "ASTER-wh"
namemap["wastrid_s"] = "ASTRID-ws"
namemap["wtreeqmc_wf_n2"] = "TQMC-n2"
namemap["wtreeqmc_wn_n2"] = "TQMC-wn_n2"
namemap["wtreeqmc_ws_n2"] = "TQMC-ws_n2"
namemap["wtreeqmc_wh_n2"] = "TQMC-wh_n2"
namemap["wtreeqmc_wh_n1"] = "TQMC-wh_n1"
namemap["wtreeqmc_wh_n0"] = "TQMC-wh_n0"

nbpss = [200, 400, 800, 1600]
ngens = [50, 200, 500, 1000]
supps = ["abayes", "bs"]
repls = range(1, 51)

for nbps in nbpss:
    for ngen in ngens:
        for repl in repls:
            print("%d %d %d" % (nbps, ngen, repl))
            for supp in supps:
                # Process quartet score for true species tree
                xqsc_df = qsc_df[(qsc_df["NBPS"] == nbps) &
                                 (qsc_df["NGEN"] == ngen) &
                                 (qsc_df["REPL"] == repl) &
                                 (qsc_df["SUPP"] == supp) &
                                 (qsc_df["MTHD"] == "true_stree")]

                if xqsc_df.shape[0] != 2:
                    sys.exit("  0 ERROR - true_stree!\n")

                true_qswn = xqsc_df[(xqsc_df["QSUP"] == "qsupp-wn")].QSCR.values[0]
                true_qswh = xqsc_df[(xqsc_df["QSUP"] == "qsupp-wh")].QSCR.values[0]

                for mthd in mthds:
                
                    # Process species tree error
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

                    qswn = "NA"
                    qswh = "NA"
                    if (mthd == "aster_h") or (mthd == "wtreeqmc_wh_n2"):
                        # Process quartet score
                        xqsc_df = qsc_df[(qsc_df["NBPS"] == nbps) &
                                         (qsc_df["NGEN"] == ngen) &
                                         (qsc_df["REPL"] == repl) &
                                         (qsc_df["SUPP"] == supp) &
                                         (qsc_df["MTHD"] == mthd)]

                        if xqsc_df.shape[0] != 2:
                            sys.exit("  2 ERROR - %s!\n" % mthd)

                        qswn = xqsc_df[(xqsc_df["QSUP"] == "qsupp-wn")].QSCR.values[0]
                        qswh = xqsc_df[(xqsc_df["QSUP"] == "qsupp-wh")].QSCR.values[0]

                    row = {}
                    row["NBPS"] = nbps
                    row["NGEN"] = ngen
                    row["REPL"] = repl
                    row["SUPP"] = supp
                    row["MTHD"] = namemap[mthd]
                    row["SEFN"] = sefn
                    row["SEFP"] = sefp
                    row["SERF"] = serf
                    row["QSWN"] = qswn
                    row["QSWH"] = qswh
                    row["TRUE_QSWN"] = true_qswn
                    row["TRUE_QSWH"] = true_qswh
                    rows.append(row)

df = pandas.DataFrame(rows, columns=cols)
df.to_csv("data-all-error-and-qscore.csv",
          sep=',', na_rep="NA",header=True, index=False)
