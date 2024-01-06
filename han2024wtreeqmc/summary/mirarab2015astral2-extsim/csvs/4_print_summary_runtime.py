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

mthds = ["ASTRID-ws",
         "TQMC-wh_n2",
         "ASTER-wh (1 thread)",
         "ASTER-wh (16 threads)"]

names = ["            ASTRID-ws",
         "           TQMC-wh_n2",
         "ASTER-wh   (1 thread)",
         "ASTER-wh (16 threads)"]

supps = ["sh"]
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

    df = pandas.read_csv("data-" + do + "-error-and-runtime.csv", keep_default_na=True)

    for ntax in ntaxs:
        for hght in hghts:
            for rate in rates:
                for ngen in ngens:
                    sys.stdout.write("Model : %s ntax, %s ILS, %s speciation, %d genes\n" % (ntax, hght, rate, ngen))
                    xdf = df[(df["NTAX"] == ntax) &
                             (df["ILSL"] == hght) &
                             (df["SPEC"] == rate) &
                             (df["NGEN"] == ngen) ]

                    for name, mthd in zip(names, mthds):
                        sys.stdout.write("  %s : " % name)

                        nrepl = len(xdf[(xdf["MTHD"] == mthds[0]) &
                                        (xdf["SUPP"] == supps[0])].REPL.values)

                        for supp in supps:
                            ydf = xdf[(xdf["MTHD"] == mthd) & (xdf["SUPP"] == supp)]

                            data = ydf.SECS.values
                            if nrepl != len(data) or nrepl > 50:
                                   sys.exit("wrong number of replicates")
                            fnravg = numpy.nanmean(data)
                            fnravg = numpy.round(numpy.round(fnravg, 5), 4)

                            sys.stdout.write("%1.4f (%s) " % (fnravg, supp))

                        sys.stdout.write("\n")
                    sys.stdout.write("\n")

"""
Increasing number of taxa
Model : 10 ntax, medium ILS, shallow speciation, 1000 genes
              ASTRID-ws : 0.0520 (sh) 
             TQMC-wh_n2 : 0.3311 (sh) 
  ASTER-wh   (1 thread) : 0.6520 (sh) 
  ASTER-wh (16 threads) : 0.1241 (sh) 

Model : 50 ntax, medium ILS, shallow speciation, 1000 genes
              ASTRID-ws : 0.1514 (sh) 
             TQMC-wh_n2 : 7.3293 (sh) 
  ASTER-wh   (1 thread) : 31.1533 (sh) 
  ASTER-wh (16 threads) : 2.4532 (sh) 

Model : 100 ntax, medium ILS, shallow speciation, 1000 genes
              ASTRID-ws : 0.1328 (sh) 
             TQMC-wh_n2 : 39.6528 (sh) 
  ASTER-wh   (1 thread) : 126.5367 (sh) 
  ASTER-wh (16 threads) : 9.3997 (sh) 

Model : 200 ntax, medium ILS, shallow speciation, 1000 genes
              ASTRID-ws : 0.2691 (sh) 
             TQMC-wh_n2 : 162.9790 (sh) 
  ASTER-wh   (1 thread) : 644.9468 (sh) 
  ASTER-wh (16 threads) : 36.8186 (sh) 

Model : 500 ntax, medium ILS, shallow speciation, 1000 genes
              ASTRID-ws : 1.4600 (sh) 
             TQMC-wh_n2 : 1059.7883 (sh) 
  ASTER-wh   (1 thread) : 4923.5136 (sh) 
  ASTER-wh (16 threads) : 262.0027 (sh) 

Model : 1000 ntax, medium ILS, shallow speciation, 1000 genes
              ASTRID-ws : 7.0398 (sh) 
             TQMC-wh_n2 : 4321.9724 (sh) 
  ASTER-wh   (1 thread) : 21975.3468 (sh) 
  ASTER-wh (16 threads) : 1077.5826 (sh) 

Increasing ILS
Model : 200 ntax, low ILS, deep speciation, 1000 genes
              ASTRID-ws : 0.2753 (sh) 
             TQMC-wh_n2 : 166.6232 (sh) 
  ASTER-wh   (1 thread) : 605.6976 (sh) 
  ASTER-wh (16 threads) : 33.3675 (sh) 

Model : 200 ntax, low ILS, shallow speciation, 1000 genes
              ASTRID-ws : 0.2671 (sh) 
             TQMC-wh_n2 : 157.2384 (sh) 
  ASTER-wh   (1 thread) : 624.2139 (sh) 
  ASTER-wh (16 threads) : 40.5089 (sh) 

Model : 200 ntax, medium ILS, deep speciation, 1000 genes
              ASTRID-ws : 0.2747 (sh) 
             TQMC-wh_n2 : 177.3989 (sh) 
  ASTER-wh   (1 thread) : 652.3588 (sh) 
  ASTER-wh (16 threads) : 37.1433 (sh) 

Model : 200 ntax, medium ILS, shallow speciation, 1000 genes
              ASTRID-ws : 0.2691 (sh) 
             TQMC-wh_n2 : 162.9790 (sh) 
  ASTER-wh   (1 thread) : 644.9468 (sh) 
  ASTER-wh (16 threads) : 36.8186 (sh) 

Model : 200 ntax, high ILS, deep speciation, 1000 genes
              ASTRID-ws : 0.2899 (sh) 
             TQMC-wh_n2 : 191.3823 (sh) 
  ASTER-wh   (1 thread) : 915.9786 (sh) 
  ASTER-wh (16 threads) : 47.6838 (sh) 

Model : 200 ntax, high ILS, shallow speciation, 1000 genes
              ASTRID-ws : 0.2850 (sh) 
             TQMC-wh_n2 : 185.9239 (sh) 
  ASTER-wh   (1 thread) : 850.9421 (sh) 
  ASTER-wh (16 threads) : 58.5240 (sh) 
"""