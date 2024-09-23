import csv
import pandas
import numpy
import sys

"""
Previous studies removed:
+ 3 remove 1 replicate (i.e. 41) from 10-taxon
+ 6 remove 2 replicates (i.e. 21, 41) from 50-taxon (note 27 was close with 590ish genes)
+ 6 remove 2 replicates (i.e. 8 and 47) from 100-taxon
+ 9 remove 3 replicates (i.e. 8, 15, 49) from 200-taxon/500K/1e-06

Concatenation (CA-ML) tree was also missing on 
+ 3 remove 1 replicate (i.e. 27) for 50-taxon, related to above
+ 1 remove 1 replicate (i.e. 12) for 10-taxon with 1000 genes
+ Also CA-ML wasn't run on any 1000-taxon with 1000 genes

TREE-QMC returned polytomies on 
+ replicate 41 for 10-taxon with 50 genes (but already removed)
+ replicate 8 for 200-taxon/500K/1e-06 with 200 genes (but already removed)

Failures due to wall clock time for single-threaded ASTER-h:
+ 2 remove 2 replicates (8 and 38) from 1000-taxon with 1000 genes
+ 1 remove 1 replicate (21) from 200-taxon/500K/1e-07 with 1000 genes
"""

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


namemap = {}
namemap["caml"] = "CA-ML"
namemap["asteroid"] = "ASTEROID"
namemap["aster_h_t16"] = "ASTER-wh"
namemap["wastrid_s"] = "ASTRID-ws"
namemap["wtreeqmc_wf_n2"] = "TQMC-n2"
namemap["wtreeqmc_wh_n2"] = "TQMC-wh_n2"

ste_df = pandas.read_csv("all_species_tree_error.csv.gz", keep_default_na=False, compression="gzip")
qsc_df = pandas.read_csv("all_quartet_score.csv.gz", keep_default_na=False, compression="gzip")
brs_df = pandas.read_csv("all_branch_support.csv.gz", keep_default_na=False, compression="gzip")

cols = ["NTAX", "ILSL", "SPEC", "REPL", "SUPP", "NGEN", 
        "MTHD",
        "SEFN", "SEFNR",
        "QSCR",
        "AVG_LPP", "MED_LPP",
        "AVG_QQS", "MED_QQS"]

ngens = [50, 200, 1000]

for do in ["varyils", "varyntax"]:
    if do == "varyntax":
        sys.stdout.write("Increasing number of taxa\n")
        ntaxs    = [10, 50, 100, 200, 500, 1000]
        hghts    = [2000000]
        hghtlabs = ["medium"]
        rates    = [0.000001]
        ratelabs = ["shallow"]

    elif do == "varyils":
        sys.stdout.write("Increasing ILS\n")
        ntaxs    = [200]
        hghts    = [10000000, 2000000, 500000]
        hghtlabs = ["low", "medium", "high"]
        rates    = [0.0000001, 0.000001]
        ratelabs = ["deep", "shallow"]

    rows = []

    for ntax in ntaxs:
        for hghtlab, hght in zip(hghtlabs, hghts):
            for ratelab, rate in zip(ratelabs, rates):
                for ngen in ngens:
                    # Pick replicates
                    repls = [repl for repl in range(1, 51)]
                    if ntax == 10:
                        repls.remove(41)
                        if ngen == 1000:
                            repls.remove(12)  # additional because of concatenation
                    elif ntax == 50:
                        repls.remove(21)
                        repls.remove(27)      # additional because of concatenation
                        repls.remove(41)
                    elif ntax == 100:
                        repls.remove(8)
                        repls.remove(47)
                    elif (ntax == 1000) and (ngen == 1000):
                            repls.remove(8)  # additional because of aster single threaded failure
                            repls.remove(38) # additional because of aster single threaded failure
                    elif (hght == 500000):
                        if (rate == 1e-6):
                            repls.remove(8)
                            repls.remove(15)
                            repls.remove(49)
                        else:
                            if ngen == 1000:
                                repls.remove(21) # additional because of aster single threaded failure

                    # Pick methods
                    mthds = ["caml", 
                             "asteroid",
                             "wastrid_s",
                             "wtreeqmc_wf_n2",
                             "wtreeqmc_wh_n2",
                             "aster_h_t16"]

                    for repl in repls:
                        print("%d taxa, %d %g, %d genes : repl %d" % (ntax, hght, rate, ngen, repl))
                        for mthd in mthds:
                            if mthd == "caml":
                                supps = ["sh"]
                            else:
                                supps = ["abayes"]

                            for supp in supps:
                                # Process species tree error
                                sefn = "NA"
                                sefnr = "NA"
                                if not ((mthd == "caml") and (ntax == 1000) and (ngen == 1000)):
                                    xste_df = ste_df[(ste_df["NTAX"] == ntax) &
                                                     (ste_df["HGHT"] == hght) &
                                                     (ste_df["RATE"] == rate) &
                                                     (ste_df["REPL"] == repl) &
                                                     (ste_df["NGEN"] == ngen) &
                                                     (ste_df["MTHD"] == mthd) &
                                                     (ste_df["SUPP"] == supp)]

                                    if xste_df.shape[0] != 1:
                                        sys.exit("  PROBLEM 1 - %s %s" % (mthd, supp))
                                    else:
                                        int1 = int(xste_df.I1.values[0])
                                        int2 = int(xste_df.I2.values[0])
                                        if int1 != int2:
                                            print("  POLY - %s %s!" % (mthd, supp))
                                        sefn = int(xste_df.FN.values[0])
                                        sefnr = float(sefn) / int1

                                # Process quartet score
                                qscr = "NA"
                                if (mthd == "aster_h_t16") or (mthd == "wtreeqmc_wh_n2"):
                                    xqsc_df = qsc_df[(qsc_df["NTAX"] == ntax) &
                                                     (qsc_df["HGHT"] == hght) &
                                                     (qsc_df["RATE"] == rate) &
                                                     (qsc_df["REPL"] == repl) &
                                                     (qsc_df["NGEN"] == ngen) &
                                                     (qsc_df["MTHD"] == mthd) &
                                                     (qsc_df["SUPP"] == supp)]

                                    if xqsc_df.shape[0] != 1:
                                        sys.exit("  PROBLEM 2 - %s %s" % (mthd, supp))
                                    else:
                                        qscr = float(xqsc_df.QS.values[0])

                                # Process branch support
                                avg_lpp = "NA"
                                med_lpp = "NA"
                                avg_qqs = "NA"
                                med_qqs = "NA"
                                if (mthd == "aster_h_t16") or (mthd == "wtreeqmc_wh_n2"):
                                    xbrs_df = brs_df[(brs_df["NTAX"] == ntax) &
                                                     (brs_df["HGHT"] == hght) &
                                                     (brs_df["RATE"] == rate) &
                                                     (brs_df["REPL"] == repl) &
                                                     (brs_df["NGEN"] == ngen) &
                                                     (brs_df["SUPP"] == supp)]

                                    if xbrs_df.shape[0] != 1:
                                        sys.exit("  3 BRANCH SUPPORT - %s!\n" % mthd)

                                    if mthd == "aster_h_t16":
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
                                row["NTAX"] = ntax
                                row["ILSL"] = hghtlab
                                row["SPEC"] = ratelab
                                row["REPL"] = repl
                                row["SUPP"] = supp
                                row["NGEN"] = ngen
                                row["MTHD"] = namemap[mthd]
                                row["SEFN"] = sefn
                                row["SEFNR"] = sefnr
                                row["QSCR"] = qscr
                                row["AVG_LPP"] = avg_lpp
                                row["MED_LPP"] = med_lpp
                                row["AVG_QQS"] = avg_qqs
                                row["MED_QQS"] = med_qqs
                                rows.append(row)

    ydf = pandas.DataFrame(rows, columns=cols)
    ydf.to_csv("data-" + do + "-error-and-qscore.csv",
               sep=',', na_rep="NA",header=True, index=False,
               quoting=csv.QUOTE_NONNUMERIC)
