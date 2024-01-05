import pandas
import numpy
import sys

### DATA FROM FIGURE 2B in WASTRAL paper (varying ILS level and speciation rate)
### ALSO SEE FIGURE 3 in WASTRID paper
### DATA for varying number of taxa is also included

### IMPORTANT: We ran Asteroid and found that it had 100% error!
###            We looked at a few cases and found that the error was much better when using
###            -n option, which turns off the missing data correction. However, this is the
###            same as running ASTRID, and wASTRID outperforms ASTRID in these data.
###            Therefore, we don't include the results Asteroid results going forward.

mthds = ["treeqmc_n2_v1.0.0",
         "wtreeqmc_wn_n2",
         "wtreeqmc_ws_n2",
         "wtreeqmc_wh_n2"]

names = ["treeqmc_n2_v1.0.0",
         "   wtreeqmc_wn_n2",
         "   wtreeqmc_ws_n2",
         "   wtreeqmc_wh_n2"]

ngens = [1000]

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

    df = pandas.read_csv("data-" + do + "-error-vs-original.csv", keep_default_na=True)

    for ntax in ntaxs:
        for hght in hghts:
            for rate in rates:
                for ngen in ngens:
                    sys.stdout.write("Model : %s ntax, %s hght, %s rate, %d genes\n" % (ntax, hght, rate, ngen))
                    xdf = df[(df["NTAX"] == ntax) &
                             (df["HGHT"] == hght) &
                             (df["RATE"] == rate) &
                             (df["NGEN"] == ngen)]

                    for mthd, name in zip(mthds, names):
                        if (mthd == "treeqmc_n2_v1.0.0") or (mthd == "wtreeqmc_wn_n2"):
                            supps = ["none"]
                            shift = "            "
                        else:
                            supps = ["sh", "abayes"]
                            shift = ""

                        sys.stdout.write("  %s : " % name)

                        for supp in supps:
                            ydf = xdf[(xdf["MTHD"] == mthd)  & (xdf["SUPP"] == supp)]

                            data = ydf.SEFNR.values
                            fnravg = numpy.nanmean(data)
                            fnravg = numpy.round(numpy.round(fnravg, 5), 4)

                            sys.stdout.write("%s%1.4f (%s) " % (shift, fnravg, supp))

                        sys.stdout.write("\n")
                    sys.stdout.write("\n")


"""
Increasing number of taxa
Model : 10 ntax, 2000000 hght, 1e-06 rate, 1000 genes
  treeqmc_n2_v1.0.0 :             0.0156 (none) 
     wtreeqmc_wn_n2 :             0.0156 (none) 
     wtreeqmc_ws_n2 : 0.0182 (sh) 0.0156 (abayes) 
     wtreeqmc_wh_n2 : 0.0182 (sh) 0.0156 (abayes) 

Model : 50 ntax, 2000000 hght, 1e-06 rate, 1000 genes
  treeqmc_n2_v1.0.0 :             0.0248 (none) 
     wtreeqmc_wn_n2 :             0.0253 (none) 
     wtreeqmc_ws_n2 : 0.0208 (sh) 0.0177 (abayes) 
     wtreeqmc_wh_n2 : 0.0226 (sh) 0.0177 (abayes) 

Model : 100 ntax, 2000000 hght, 1e-06 rate, 1000 genes
  treeqmc_n2_v1.0.0 :             0.0306 (none) 
     wtreeqmc_wn_n2 :             0.0304 (none) 
     wtreeqmc_ws_n2 : 0.0276 (sh) 0.0221 (abayes) 
     wtreeqmc_wh_n2 : 0.0255 (sh) 0.0198 (abayes) 

Model : 200 ntax, 2000000 hght, 1e-06 rate, 1000 genes
  treeqmc_n2_v1.0.0 :             0.0315 (none) 
     wtreeqmc_wn_n2 :             0.0310 (none) 
     wtreeqmc_ws_n2 : 0.0279 (sh) 0.0236 (abayes) 
     wtreeqmc_wh_n2 : 0.0262 (sh) 0.0227 (abayes) 

Model : 500 ntax, 2000000 hght, 1e-06 rate, 1000 genes
  treeqmc_n2_v1.0.0 :             0.0280 (none) 
     wtreeqmc_wn_n2 :             0.0276 (none) 
     wtreeqmc_ws_n2 : 0.0258 (sh) 0.0237 (abayes) 
     wtreeqmc_wh_n2 : 0.0227 (sh) 0.0203 (abayes) 

Model : 1000 ntax, 2000000 hght, 1e-06 rate, 1000 genes
  treeqmc_n2_v1.0.0 :             0.0416 (none) 
     wtreeqmc_wn_n2 :             0.0411 (none) 
     wtreeqmc_ws_n2 : 0.0369 (sh) 0.0318 (abayes) 
     wtreeqmc_wh_n2 : 0.0315 (sh) 0.0258 (abayes) 

Increasing ILS
Model : 200 ntax, 10000000 hght, 1e-07 rate, 1000 genes
  treeqmc_n2_v1.0.0 :             0.0450 (none) 
     wtreeqmc_wn_n2 :             0.0446 (none) 
     wtreeqmc_ws_n2 : 0.0439 (sh) 0.0401 (abayes) 
     wtreeqmc_wh_n2 : 0.0332 (sh) 0.0248 (abayes) 

Model : 200 ntax, 10000000 hght, 1e-06 rate, 1000 genes
  treeqmc_n2_v1.0.0 :             0.0183 (none) 
     wtreeqmc_wn_n2 :             0.0175 (none) 
     wtreeqmc_ws_n2 : 0.0157 (sh) 0.0120 (abayes) 
     wtreeqmc_wh_n2 : 0.0148 (sh) 0.0112 (abayes) 

Model : 200 ntax, 2000000 hght, 1e-07 rate, 1000 genes
  treeqmc_n2_v1.0.0 :             0.0396 (none) 
     wtreeqmc_wn_n2 :             0.0397 (none) 
     wtreeqmc_ws_n2 : 0.0377 (sh) 0.0357 (abayes) 
     wtreeqmc_wh_n2 : 0.0341 (sh) 0.0309 (abayes) 

Model : 200 ntax, 2000000 hght, 1e-06 rate, 1000 genes
  treeqmc_n2_v1.0.0 :             0.0315 (none) 
     wtreeqmc_wn_n2 :             0.0310 (none) 
     wtreeqmc_ws_n2 : 0.0279 (sh) 0.0236 (abayes) 
     wtreeqmc_wh_n2 : 0.0262 (sh) 0.0227 (abayes) 

Model : 200 ntax, 500000 hght, 1e-07 rate, 1000 genes
  treeqmc_n2_v1.0.0 :             0.0628 (none) 
     wtreeqmc_wn_n2 :             0.0609 (none) 
     wtreeqmc_ws_n2 : 0.0562 (sh) 0.0492 (abayes) 
     wtreeqmc_wh_n2 : 0.0554 (sh) 0.0487 (abayes) 

Model : 200 ntax, 500000 hght, 1e-06 rate, 1000 genes
  treeqmc_n2_v1.0.0 :             0.0498 (none) 
     wtreeqmc_wn_n2 :             0.0493 (none) 
     wtreeqmc_ws_n2 : 0.0459 (sh) 0.0404 (abayes) 
     wtreeqmc_wh_n2 : 0.0436 (sh) 0.0377 (abayes) 
"""