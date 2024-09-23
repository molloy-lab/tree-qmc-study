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

mthds = ["ASTER-wh (1 thread)",
         "TQMC-wh_n2",
         "ASTER-wh (16 threads)",
         "TQMC-n2",
         "ASTRID-ws"]

namemap = {}
namemap["ASTRID-ws"]             = "            ASTRID-ws"
namemap["ASTER-wh (1 thread)"]   = "  ASTER-wh (1 thread)"
namemap["ASTER-wh (16 threads)"] = "ASTER-wh (16 threads)"
namemap["TQMC-wh_n2"]            = "           TQMC-wh_n2"
namemap["TQMC-n2"]               = "              TQMC-n2"

ngens = [1000]

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

    df = pandas.read_csv("data-" + do + "-runtime.csv", keep_default_na=True)

    for ntax in ntaxs:
        for hght in hghts:
            for rate in rates:
                for ngen in ngens:
                    sys.stdout.write("Model : %s ntax, %s ILS, %s speciation, %d genes\n" % (ntax, hght, rate, ngen))
                    xdf = df[(df["NTAX"] == ntax) &
                             (df["ILSL"] == hght) &
                             (df["SPEC"] == rate) &
                             (df["NGEN"] == ngen) ]

                    nrepl = len(xdf[(xdf["MTHD"] == mthds[0]) & (xdf["SUPP"] == "abayes")].REPL.values)

                    for mthd in mthds:
                        ydf = xdf[(xdf["MTHD"] == mthd) & (xdf["SUPP"] == "abayes")]

                        data = ydf.SECS.values / (60 * 60)
                        if nrepl != len(data) or nrepl > 50:
                            sys.exit("wrong number of replicates")
                        fnravg = numpy.nanmean(data)
                        fnravg = numpy.round(numpy.round(fnravg, 5), 4)

                        sys.stdout.write("%s : %1.4f hours\n" % (namemap[mthd], fnravg))

                    sys.stdout.write("\n")
                sys.stdout.write("\n")

"""
Increasing number of taxa
Model : 10 ntax, medium ILS, shallow speciation, 1000 genes
  ASTER-wh (1 thread) : 0.0002 hours
           TQMC-wh_n2 : 0.0001 hours
ASTER-wh (16 threads) : 0.0000 hours
              TQMC-n2 : 0.0001 hours
            ASTRID-ws : 0.0000 hours


Model : 50 ntax, medium ILS, shallow speciation, 1000 genes
  ASTER-wh (1 thread) : 0.0096 hours
           TQMC-wh_n2 : 0.0021 hours
ASTER-wh (16 threads) : 0.0007 hours
              TQMC-n2 : 0.0009 hours
            ASTRID-ws : 0.0000 hours


Model : 100 ntax, medium ILS, shallow speciation, 1000 genes
  ASTER-wh (1 thread) : 0.0366 hours
           TQMC-wh_n2 : 0.0112 hours
ASTER-wh (16 threads) : 0.0030 hours
              TQMC-n2 : 0.0038 hours
            ASTRID-ws : 0.0000 hours


Model : 200 ntax, medium ILS, shallow speciation, 1000 genes
  ASTER-wh (1 thread) : 0.1836 hours
           TQMC-wh_n2 : 0.0464 hours
ASTER-wh (16 threads) : 0.0114 hours
              TQMC-n2 : 0.0144 hours
            ASTRID-ws : 0.0001 hours


Model : 500 ntax, medium ILS, shallow speciation, 1000 genes
  ASTER-wh (1 thread) : 1.4294 hours
           TQMC-wh_n2 : 0.3011 hours
ASTER-wh (16 threads) : 0.0758 hours
              TQMC-n2 : 0.0896 hours
            ASTRID-ws : 0.0005 hours


Model : 1000 ntax, medium ILS, shallow speciation, 1000 genes
  ASTER-wh (1 thread) : 7.9886 hours
           TQMC-wh_n2 : 1.2416 hours
ASTER-wh (16 threads) : 0.3535 hours
              TQMC-n2 : 0.3620 hours
            ASTRID-ws : 0.0020 hours


Increasing ILS
Model : 200 ntax, low ILS, deep speciation, 1000 genes
  ASTER-wh (1 thread) : 0.1945 hours
           TQMC-wh_n2 : 0.0474 hours
ASTER-wh (16 threads) : 0.0115 hours
              TQMC-n2 : 0.0144 hours
            ASTRID-ws : 0.0001 hours


Model : 200 ntax, low ILS, shallow speciation, 1000 genes
  ASTER-wh (1 thread) : 0.1793 hours
           TQMC-wh_n2 : 0.0444 hours
ASTER-wh (16 threads) : 0.0102 hours
              TQMC-n2 : 0.0140 hours
            ASTRID-ws : 0.0001 hours


Model : 200 ntax, medium ILS, deep speciation, 1000 genes
  ASTER-wh (1 thread) : 0.2040 hours
           TQMC-wh_n2 : 0.0504 hours
ASTER-wh (16 threads) : 0.0122 hours
              TQMC-n2 : 0.0155 hours
            ASTRID-ws : 0.0001 hours


Model : 200 ntax, medium ILS, shallow speciation, 1000 genes
  ASTER-wh (1 thread) : 0.1836 hours
           TQMC-wh_n2 : 0.0464 hours
ASTER-wh (16 threads) : 0.0114 hours
              TQMC-n2 : 0.0144 hours
            ASTRID-ws : 0.0001 hours


Model : 200 ntax, high ILS, deep speciation, 1000 genes
  ASTER-wh (1 thread) : 0.2653 hours
           TQMC-wh_n2 : 0.0560 hours
ASTER-wh (16 threads) : 0.0150 hours
              TQMC-n2 : 0.0170 hours
            ASTRID-ws : 0.0001 hours


Model : 200 ntax, high ILS, shallow speciation, 1000 genes
  ASTER-wh (1 thread) : 0.2415 hours
           TQMC-wh_n2 : 0.0535 hours
ASTER-wh (16 threads) : 0.0150 hours
              TQMC-n2 : 0.0163 hours
            ASTRID-ws : 0.0001 hours
"""