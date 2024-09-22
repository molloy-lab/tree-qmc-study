import pandas
import numpy
import sys

sys.exit("DONE RUNNING")

def string_to_float_array(fptext, tptext):
    ps = []

    if fptext == '':
        pass
    else:
        fps = fptext.split(',')
        ps += fps

    if tptext == '':
        pass
    else:
        tps = tptext.split(',')
        ps += tps

    keep = []
    for x in ps:
        if x != "NA":
            y = float(x)

            # Some weird outputs can occur when EN is small
            if (y < 0) or (y > 1):
                y = 0.0

            keep.append(float(y))

    return keep


ste_df = pandas.read_csv("all_species_tree_error.csv.gz", keep_default_na=False, compression="gzip")
qsc_df = pandas.read_csv("all_quartet_score.csv.gz", keep_default_na=False, compression="gzip")
brs_df = pandas.read_csv("all_branch_support.csv.gz", keep_default_na=False, compression="gzip")

cols = ["NBPS", "NGEN", "REPL", "SUPP",
        "MTHD", 
        "SEFN", "SEFP", "SERF",
        "QSCR",
        "AVG_LPP", "MED_LPP",
        "AVG_QQS", "MED_QQS"]
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
                #xqsc_df = qsc_df[(qsc_df["NBPS"] == nbps) &
                #                 (qsc_df["NGEN"] == ngen) &
                #                 (qsc_df["REPL"] == repl) &
                #                 (qsc_df["SUPP"] == supp) &
                #                 (qsc_df["MTHD"] == "true_stree")]

                #if xqsc_df.shape[0] != 2:
                #    sys.exit("  0 ERROR - true_stree!\n")

                #true_qswn = xqsc_df[(xqsc_df["QSUP"] == "qsupp-wn")].QSCR.values[0]
                #true_qswh = xqsc_df[(xqsc_df["QSUP"] == "qsupp-wh")].QSCR.values[0]

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

                    # Process quartet score
                    qscr = "NA"
                    if (mthd == "aster_h") or (mthd == "wtreeqmc_wh_n2"):
                        xqsc_df = qsc_df[(qsc_df["NBPS"] == nbps) &
                                         (qsc_df["NGEN"] == ngen) &
                                         (qsc_df["REPL"] == repl) &
                                         (qsc_df["SUPP"] == supp) &
                                         (qsc_df["MTHD"] == mthd)]

                        if xqsc_df.shape[0] != 1:
                            sys.exit("  2 ERROR - %s!\n" % mthd)

                        qscr = float(xqsc_df.QS.values[0])

                    # Process average and median branch support
                    avg_lpp = "NA"
                    med_lpp = "NA"
                    avg_qqs = "NA"
                    med_qqs = "NA"
                    if (mthd == "aster_h") or (mthd == "wtreeqmc_wh_n2"):
                        xbrs_df = brs_df[(brs_df["NBPS"] == nbps) &
                                         (brs_df["NGEN"] == ngen) &
                                         (brs_df["REPL"] == repl) &
                                         (brs_df["SUPP"] == supp)]

                        if xbrs_df.shape[0] != 1:
                            sys.exit("  3 BRANCH SUPPORT - %s!\n" % mthd)

                        if mthd == "aster_h":
                            # estimated tree 1
                            lpp = string_to_float_array(xbrs_df["FP_E1_PP"].values[0], xbrs_df["TP_E1_PP"].values[0])
                            qqs = string_to_float_array(xbrs_df["FP_E1_QS"].values[0], xbrs_df["TP_E1_QS"].values[0])
                        elif mthd == "wtreeqmc_wh_n2":
                            # estimated tree 2
                            lpp = string_to_float_array(xbrs_df["FP_E2_PP"].values[0], xbrs_df["TP_E2_PP"].values[0])
                            qqs = string_to_float_array(xbrs_df["FP_E2_QS"].values[0], xbrs_df["TP_E2_QS"].values[0])
                        else:
                            sys.exit("Branch support was not computed and compared for other methods!")
                        
                        avg_lpp = numpy.mean(lpp)
                        med_lpp = numpy.median(lpp)
                        
                        avg_qqs = numpy.mean(qqs)
                        med_qqs = numpy.median(qqs)

                    row = {}
                    row["NBPS"] = nbps
                    row["NGEN"] = ngen
                    row["REPL"] = repl
                    row["SUPP"] = supp
                    row["MTHD"] = namemap[mthd]
                    row["SEFN"] = sefn
                    row["SEFP"] = sefp
                    row["SERF"] = serf
                    row["QSCR"] = qscr
                    row["AVG_LPP"] = avg_lpp
                    row["MED_LPP"] = med_lpp
                    row["AVG_QQS"] = avg_qqs
                    row["MED_QQS"] = med_qqs
                    rows.append(row)

df = pandas.DataFrame(rows, columns=cols)
df.to_csv("data-all-error-and-qscore.csv",
          sep=',', na_rep="NA",header=True, index=False)
