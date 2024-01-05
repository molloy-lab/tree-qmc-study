import pandas
import numpy
import sys

"""
Previous studies removed:
+ remove 1 replicate (i.e. 41) from 10-taxon
+ remove 3 replicates (i.e. 21 and 41) from 50-taxon (next is 27 with 590ish genes)
+ remove 2 replicates (i.e. 8 and 47) from 100-taxon
"""

#sys.exit("DONE RUNNING")

df1 = pandas.read_csv("all_species_tree_error.csv",
                      keep_default_na=False)
df2 = pandas.read_csv("caml_species_tree_error.csv",
                      keep_default_na=False)

ste_df = pandas.concat([df1, df2])

cols = ["NTAX", "HGHT", "RATE", "REPL",
        "SUPP", "NGEN", "MTHD", "SEFN", "SEFNR"]

ngens = [50, 200, 1000]
mthds = ["wastrid_s", "asteroid", "aster_h_t16", "wtreeqmc_wh_n2",
         "wtreeqmc_wh_n1", "wtreeqmc_wh_n0", "wtreeqmc_ws_n2", "wtreeqmc_wn_n2",
         "caml"]

for do in ["varyntax", "varyils"]: 
    if do == "varyntax":
        sys.stdout.write("Increasing number of taxa\n")
        ntaxs = [10, 50, 100, 200, 500, 1000]
        hghts = [2000000]
        rates = [0.000001]
        
    elif do == "varyils":
        sys.stdout.write("Increasing ILS\n")
        ntaxs = [200]
        hghts = [10000000, 2000000, 500000]
        rates = [0.0000001, 0.000001]

    rows = []

    for ntax in ntaxs:
        for hght in hghts:
            for rate in rates:
                for ngen in ngens:
                    repls = [repl for repl in range(1, 51)]
                    if ntax == 10:
                        repls.remove(41)
                        if ngen == 1000:
                            repls.remove(12)  # additional
                    elif ntax == 50:
                        repls.remove(21)
                        repls.remove(27)  # additional
                        repls.remove(41)
                    elif ntax == 100:
                        repls.remove(8)
                        repls.remove(47)
                    elif (hght == 500000) and (rate == 1e-6):
                        repls.remove(8)
                        repls.remove(15)
                        repls.remove(49)
                    for repl in repls:
                        print("%d taxa, %d %g, %d genes : repl %d" % (ntax, hght, rate, ngen, repl))
                        for mthd in mthds:
                            if (mthd == "asteroid") or (mthd == "wtreeqmc_wn_n2") or (mthd == "caml"):
                                supps = ["sh"]
                            else:
                                supps = ["sh", "abayes"]

                            for supp in supps:
                                xste_df = ste_df[(ste_df["NTAX"] == ntax) &
                                         (ste_df["HGHT"] == hght) &
                                         (ste_df["RATE"] == rate) &
                                         (ste_df["REPL"] == repl) & 
                                         (ste_df["SUPP"] == supp) &
                                         (ste_df["NGEN"] == ngen) &
                                         (ste_df["MTHD"] == mthd)]

                                if xste_df.shape[0] != 1:
                                    print("  PROBLEM - %s %s" % (mthd, supp))
                                    sefn = "NA"
                                    sefnr = "NA"
                                else:
                                    int1 = int(xste_df.I1.values[0])
                                    int2 = int(xste_df.I2.values[0])
                                    if int1 != int2:
                                        print("  POLY - %s %s!" % (mthd, supp))
                                    sefn = int(xste_df.FN.values[0])
                                    sefnr = float(sefn) / int1

                                row = {}
                                row["NTAX"] = ntax
                                row["HGHT"] = hght
                                row["RATE"] = rate
                                row["REPL"] = repl
                                if (mthd == "asteroid") or (mthd == "wtreeqmc_wn_n2") or (mthd == "caml"):
                                    row["SUPP"] = "none"
                                else:
                                    row["SUPP"] = supp
                                row["NGEN"] = ngen
                                row["MTHD"] = mthd
                                row["SEFN"] = sefn
                                row["SEFNR"] = sefnr
                                rows.append(row)

    df = pandas.DataFrame(rows, columns=cols)
    df.to_csv("data-" + do + "-error.csv",
              sep=',', na_rep="NA",header=True, index=False)
