import pandas
import numpy
import sys

"""
Previous studies removed:
+ remove 1 replicate (i.e. 41) from 10-taxon
+ remove 3 replicates (i.e. 21 and 41) from 50-taxon (next is 27 with 590ish genes)
+ remove 2 replicates (i.e. 8 and 47) from 100-taxon
+ remove 3 replicates (i.e. 8, 15, 49) from 200 taxon/500 K/1e-06
"""

sys.exit("DONE RUNNING")

namemap = {}
namemap["caml"] = "CA-ML"
namemap["treeqmc_n2_v1.0.0"] = "TQMC-n2-origstudy"
namemap["aster_h_t16"] = "ASTER-wh"
namemap["wastrid_s"] = "ASTRID-ws"
namemap["wtreeqmc_wn_n2"] = "TQMC-n2"
namemap["wtreeqmc_ws_n2"] = "TQMC-ws_n2"
namemap["wtreeqmc_wh_n2"] = "TQMC-wh_n2"
namemap["wtreeqmc_wh_n1"] = "TQMC-wh_n1"
namemap["wtreeqmc_wh_n0"] = "TQMC-wh_n0"

df1 = pandas.read_csv("all_species_tree_error.csv",
                      keep_default_na=False)
df2 = pandas.read_csv("caml_species_tree_error.csv",
                      keep_default_na=False)

df3 = pandas.read_csv("../../../../han2023improving/summary/mirarab2015astral2/csvs/all_species_tree_error.csv",
                      keep_default_na=False)
df3 = df3[(df3["NGEN"] == 1000) &
          (df3["MTHD"] == "treeqmc_n2_v1.0.0") &
          (df3["GTRE"] == "estimatedgenetre")]  # Get original results on 1000 gene trees
df3 = df3.rename(columns={"STRHT": "HGHT",
                          "SRATE": "RATE",
                          "GTRE": "SUPP"})
df = pandas.concat([df1, df2, df3])

cols = ["NTAX", "ILSL", "SPEC", "REPL",
        "SUPP", "NGEN", "MTHD", "SEFN", "SEFNR"]

ngens = [50, 200, 1000]

for do in ["varyils", "varyntax"]:
    if do == "varyntax":
        sys.stdout.write("Increasing number of taxa\n")
        ntaxs = [10, 50, 100, 200, 500, 1000]
        hghts = [2000000]
        hghtlabs = ["medium"]
        rates = [0.000001]
        ratelabs = ["shallow"]
        
    elif do == "varyils":
        sys.stdout.write("Increasing ILS\n")
        ntaxs = [200]
        hghts = [10000000, 2000000, 500000]
        hghtlabs = ["low", "medium", "high"]
        rates = [0.0000001, 0.000001]
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
                        repls.remove(27)  # additional because of concatenation
                        repls.remove(41)
                    elif ntax == 100:
                        repls.remove(8)
                        repls.remove(47)
                    elif (hght == 500000) and (rate == 1e-6):
                        repls.remove(8)
                        repls.remove(15)
                        repls.remove(49)

                    # Pick methods
                    mthds = ["caml",               # from prior study
                             "treeqmc_n2_v1.0.0",  # from prior study
                             "wtreeqmc_wn_n2",
                             "wtreeqmc_ws_n2",
                             "wtreeqmc_wh_n2",
                             "aster_h_t16",
                             "wastrid_s"]

                    # To suppress warnings
                    #if ngen != 1000:
                    #    mthds.remove("caml")
                    #    mthds.remove("treeqmc_n2_v1.0.0")
                    #else:
                    #    if ntax == 1000:
                    #        mthds.remove("caml")

                    for repl in repls:
                        print("%d taxa, %d %g, %d genes : repl %d" % (ntax, hght, rate, ngen, repl))
                        for mthd in mthds:
                            if (mthd == "caml") or \
                               (mthd == "treeqmc_n2_v1.0.0") or \
                               (mthd == "wtreeqmc_wn_n2"):
                                supps = ["none"]
                            else:
                                supps = ["abayes", "sh"]

                            for supp in supps:
                                xdf = df[(df["NTAX"] == ntax) &
                                         (df["HGHT"] == hght) &
                                         (df["RATE"] == rate) &
                                         (df["REPL"] == repl) &
                                         (df["NGEN"] == ngen) &
                                         (df["MTHD"] == mthd)]

                                if supp != "none":
                                    xdf = xdf[(xdf["SUPP"] == supp)]

                                if xdf.shape[0] != 1:
                                    print("  PROBLEM - %s %s" % (mthd, supp))
                                    sefn = "NA"
                                    sefnr = "NA"
                                else:
                                    int1 = int(xdf.I1.values[0])
                                    int2 = int(xdf.I2.values[0])
                                    if int1 != int2:
                                        print("  POLY - %s %s!" % (mthd, supp))
                                    sefn = int(xdf.FN.values[0])
                                    sefnr = float(sefn) / int1

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
                                rows.append(row)

    ydf = pandas.DataFrame(rows, columns=cols)
    ydf.to_csv("data-" + do + "-error.csv",
               sep=',', na_rep="NA",header=True, index=False)
