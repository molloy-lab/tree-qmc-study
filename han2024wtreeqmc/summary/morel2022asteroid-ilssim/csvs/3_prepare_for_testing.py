import pandas
import numpy
import sys

sys.exit("DONE RUNNING")

cols = ["NTAX", "NGEN", "NBPS", "BLSC", "PSIZ", "MISS", "REPL",
        "WASTRIDxFNR", "WASTRIDxFPR",
        "ASTERxFNR", "ASTERxFPR",
        "ASTEROIDxFNR", "ASTEROIDxFPR",
        "TQMCn2xFNR", "TQMCn2xFPR",
        "TQMCn2sharedxFNR", "TQMCn2sharedxFPR",
        "TQMCn1xFNR", "TQMCn1xFPR",
        "TQMCn0xFNR", "TQMCn0xFPR"]


experiments = ["varypsiz", "varyntax", "varyngen",
               "varynbps", "varyblsc", "varymiss"]

for do in experiments:
    ntaxs = [50]
    ngens = [1000]
    nbpss = [100]
    blscs = [1.0]
    psizs = [50000000]
    misss = [0.6]  # ms is always same as mf
    repls = [repl for repl in range(0, 51)]

    if do == "varypsiz":
        psizs = [10, 50000000, 100000000, 500000000, 1000000000]
    elif do == "varyntax":
        ntaxs = [25, 75, 50, 100, 125, 150]  # dup 50
    elif do == "varyngen":
        ngens = [250, 500, 1000, 2000]  # dup 1000
    elif do == "varynbps":
        nbpss = [50, 100, 200, 500]  # dup 100
    elif do == "varyblsc":
        blscs = [0.05, 0.1, 1, 10, 100, 200]  # dup 1.0
    elif do == "varymiss":
        misss = [0.5, 0.55, 0.6, 0.65, 0.7, 0.75]  # dup 0.6
    else:
        sys.exit()

    df = pandas.read_csv("data-" + do + "-error-and-timings.csv",
                         na_values='NA', keep_default_na=False)
    rows = []

    for ntax in ntaxs:
        for ngen in ngens:
            for nbps in nbpss:
                for blsc in blscs:
                    for psiz in psizs:
                        for miss in misss:
                            print("%s - %d taxa, %d genes, %d bps, %f bl scaler, %d pop size, %f missing param" \
                                  % (do, ntax, ngen, nbps, blsc, psiz, miss))

                            xdf = df[(df["NTAX"] == ntax) &
                                     (df["NGEN"] == ngen) &
                                     (df["NBPS"] == nbps) &
                                     (df["BLSC"] == blsc) &
                                     (df["PSIZ"] == psiz) & 
                                     (df["MISS"] == miss)]

                            repls = xdf[(xdf["MTHD"] == "asteroid")].REPL.values

                            for repl in repls:
                                ydf = xdf[(xdf["REPL"] == repl)]

                                #print(ydf)

                                row = {}
                                row["NTAX"] = ntax
                                row["NGEN"] = ngen
                                row["NBPS"] = nbps
                                row["BLSC"] = blsc
                                row["PSIZ"] = psiz
                                row["MISS"] = miss
                                row["REPL"] = repl

                                vals = ydf[(ydf["MTHD"] == "wastrid_vanilla")]
                                row["WASTRIDxFNR"] = vals.SEFNR.values[0]
                                row["WASTRIDxFPR"] = vals.SEFPR.values[0]

                                vals = ydf[(ydf["MTHD"] == "aster_v1.16.3.4")]
                                row["ASTERxFNR"] = vals.SEFNR.values[0]
                                row["ASTERxFPR"] = vals.SEFPR.values[0]

                                vals = ydf[(ydf["MTHD"] == "asteroid")]
                                row["ASTEROIDxFNR"] = vals.SEFNR.values[0]
                                row["ASTEROIDxFPR"] = vals.SEFPR.values[0]

                                vals = ydf[(ydf["MTHD"] == "wtreeqmc_wf_n2")]
                                row["TQMCn2xFNR"] = vals.SEFNR.values[0]
                                row["TQMCn2xFPR"] = vals.SEFPR.values[0]

                                vals = ydf[(ydf["MTHD"] == "wtreeqmc_wf_n2_shared")]
                                row["TQMCn2sharedxFNR"] = vals.SEFNR.values[0]
                                row["TQMCn2sharedxFPR"] = vals.SEFPR.values[0]

                                vals = ydf[(ydf["MTHD"] == "wtreeqmc_wf_n1")]
                                row["TQMCn1xFNR"] = vals.SEFNR.values[0]
                                row["TQMCn1xFPR"] = vals.SEFPR.values[0]

                                vals = ydf[(ydf["MTHD"] == "wtreeqmc_wf_n0")]
                                row["TQMCn0xFNR"] = vals.SEFNR.values[0]
                                row["TQMCn0xFPR"] = vals.SEFPR.values[0]

                                rows.append(row)

    df = pandas.DataFrame(rows, columns=cols)
    df.to_csv("data-" + do + "-for-testing.csv",
            sep=',', na_rep="NA",header=True, index=False)
