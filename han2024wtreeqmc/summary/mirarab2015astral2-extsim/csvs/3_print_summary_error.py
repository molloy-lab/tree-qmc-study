import pandas
import numpy
import sys
import warnings

### DATA FROM FIGURE 2B in WASTRAL paper (varying ILS level and speciation rate)
### ALSO SEE FIGURE 3 in WASTRID paper
### DATA for varying number of taxa is also included

### IMPORTANT: We ran Asteroid and found that it had 100% error!
###            We looked at a few cases and found that the error was much better when using
###            -n option, which turns off the missing data correction. However, this is the
###            same as running ASTRID, and wASTRID outperforms ASTRID in these data.
###            Therefore, we don't include the results Asteroid results going forward.

mthds = ["CA-ML",
         "ASTRID-ws",
         "ASTER-wh",
         "TQMC-wh_n2",
         "TQMC-n2",
         "ASTEROID"]

namemap = {}
namemap["CA-ML"]         = "            CA-ML"
namemap["ASTRID-ws"]     = "        ASTRID-ws"
namemap["ASTER-wh"]      = "         ASTER-wh"
namemap["TQMC-wh_n2"]    = "       TQMC-wh_n2"
namemap["TQMC-n2"]       = "          TQMC-n2"
namemap["ASTEROID"]      = "         ASTEROID"

ngens = [50, 200, 1000]

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

    df = pandas.read_csv("data-" + do + "-error-and-qscore.csv", keep_default_na=True)

    for ntax in ntaxs:
        for hght in hghts:
            for rate in rates:
                for ngen in ngens:
                    sys.stdout.write("Model : %s ntax, %s ILS, %s speciation, %d genes\n" % (ntax, hght, rate, ngen))
                    xdf = df[(df["NTAX"] == ntax) &
                             (df["ILSL"] == hght) &
                             (df["SPEC"] == rate) &
                             (df["NGEN"] == ngen)]

                    nrepl = len(xdf[(xdf["MTHD"] == "CA-ML")].REPL.values)

                    for mthd in mthds:
                        ydf = xdf[(xdf["MTHD"] == mthd)]
                        data = ydf.SEFNR.values
                        if nrepl != len(data) or nrepl > 50:
                            sys.exit("Wrong number of replicates")

                        #with warnings.catch_warnings():
                        #    warnings.simplefilter("ignore", category=RuntimeWarning)
                        fnravg = numpy.nanmean(data)
                        fnravg = numpy.round(numpy.round(fnravg, 5), 4)

                        sys.stdout.write("%s : %1.4f\n" % (namemap[mthd], fnravg))

                sys.stdout.write("\n")
            sys.stdout.write("\n")

"""
Increasing number of taxa
Model : 10 ntax, medium ILS, shallow speciation, 50 genes
            CA-ML : 0.0383
        ASTRID-ws : 0.0281
         ASTER-wh : 0.0306
       TQMC-wh_n2 : 0.0306
          TQMC-n2 : 0.0357
         ASTEROID : 0.4158
Model : 10 ntax, medium ILS, shallow speciation, 200 genes
            CA-ML : 0.0179
        ASTRID-ws : 0.0076
         ASTER-wh : 0.0076
       TQMC-wh_n2 : 0.0102
          TQMC-n2 : 0.0179
         ASTEROID : 0.4133
Model : 10 ntax, medium ILS, shallow speciation, 1000 genes
            CA-ML : 0.0208
        ASTRID-ws : 0.0156
         ASTER-wh : 0.0156
       TQMC-wh_n2 : 0.0156
          TQMC-n2 : 0.0156
         ASTEROID : 0.4141


Model : 50 ntax, medium ILS, shallow speciation, 50 genes
            CA-ML : 0.0780
        ASTRID-ws : 0.0576
         ASTER-wh : 0.0607
       TQMC-wh_n2 : 0.0585
          TQMC-n2 : 0.0709
         ASTEROID : 0.9942
Model : 50 ntax, medium ILS, shallow speciation, 200 genes
            CA-ML : 0.0452
        ASTRID-ws : 0.0355
         ASTER-wh : 0.0306
       TQMC-wh_n2 : 0.0284
          TQMC-n2 : 0.0412
         ASTEROID : 0.9938
Model : 50 ntax, medium ILS, shallow speciation, 1000 genes
            CA-ML : 0.0266
        ASTRID-ws : 0.0191
         ASTER-wh : 0.0186
       TQMC-wh_n2 : 0.0177
          TQMC-n2 : 0.0253
         ASTEROID : 0.9934


Model : 100 ntax, medium ILS, shallow speciation, 50 genes
            CA-ML : 0.0912
        ASTRID-ws : 0.0704
         ASTER-wh : 0.0716
       TQMC-wh_n2 : 0.0672
          TQMC-n2 : 0.0723
         ASTEROID : 0.9966
Model : 100 ntax, medium ILS, shallow speciation, 200 genes
            CA-ML : 0.0474
        ASTRID-ws : 0.0393
         ASTER-wh : 0.0383
       TQMC-wh_n2 : 0.0361
          TQMC-n2 : 0.0468
         ASTEROID : 0.9966
Model : 100 ntax, medium ILS, shallow speciation, 1000 genes
            CA-ML : 0.0249
        ASTRID-ws : 0.0240
         ASTER-wh : 0.0191
       TQMC-wh_n2 : 0.0198
          TQMC-n2 : 0.0304
         ASTEROID : 0.9964


Model : 200 ntax, medium ILS, shallow speciation, 50 genes
            CA-ML : 0.0922
        ASTRID-ws : 0.0728
         ASTER-wh : 0.0706
       TQMC-wh_n2 : 0.0669
          TQMC-n2 : 0.0809
         ASTEROID : 0.9985
Model : 200 ntax, medium ILS, shallow speciation, 200 genes
            CA-ML : 0.0552
        ASTRID-ws : 0.0417
         ASTER-wh : 0.0403
       TQMC-wh_n2 : 0.0382
          TQMC-n2 : 0.0475
         ASTEROID : 0.9987
Model : 200 ntax, medium ILS, shallow speciation, 1000 genes
            CA-ML : 0.0278
        ASTRID-ws : 0.0260
         ASTER-wh : 0.0241
       TQMC-wh_n2 : 0.0227
          TQMC-n2 : 0.0317
         ASTEROID : 0.9988


Model : 500 ntax, medium ILS, shallow speciation, 50 genes
            CA-ML : 0.0924
        ASTRID-ws : 0.0740
         ASTER-wh : 0.0740
       TQMC-wh_n2 : 0.0665
          TQMC-n2 : 0.0763
         ASTEROID : 0.9995
Model : 500 ntax, medium ILS, shallow speciation, 200 genes
            CA-ML : 0.0472
        ASTRID-ws : 0.0410
         ASTER-wh : 0.0400
       TQMC-wh_n2 : 0.0341
          TQMC-n2 : 0.0426
         ASTEROID : 0.9996
Model : 500 ntax, medium ILS, shallow speciation, 1000 genes
            CA-ML : 0.0233
        ASTRID-ws : 0.0261
         ASTER-wh : 0.0243
       TQMC-wh_n2 : 0.0203
          TQMC-n2 : 0.0281
         ASTEROID : 0.9996


Model : 1000 ntax, medium ILS, shallow speciation, 50 genes
            CA-ML : 0.0976
        ASTRID-ws : 0.0881
         ASTER-wh : 0.0867
       TQMC-wh_n2 : 0.0768
          TQMC-n2 : 0.1011
         ASTEROID : 0.9998
Model : 1000 ntax, medium ILS, shallow speciation, 200 genes
            CA-ML : 0.0515
        ASTRID-ws : 0.0518
         ASTER-wh : 0.0485
       TQMC-wh_n2 : 0.0418
          TQMC-n2 : 0.0588
         ASTEROID : 0.9998
Model : 1000 ntax, medium ILS, shallow speciation, 1000 genes
/Users/ekmolloy/Desktop/tree-qmc-study/han2024wtreeqmc/summary/mirarab2015astral2-extsim/csvs/3_print_summary_error.py:67: RuntimeWarning: Mean of empty slice
  fnravg = numpy.nanmean(data)
            CA-ML : nan
        ASTRID-ws : 0.0332
         ASTER-wh : 0.0305
       TQMC-wh_n2 : 0.0246
          TQMC-n2 : 0.0364
         ASTEROID : 0.9998


Increasing ILS
Model : 200 ntax, low ILS, deep speciation, 50 genes
            CA-ML : 0.0398
        ASTRID-ws : 0.0649
         ASTER-wh : 0.0583
       TQMC-wh_n2 : 0.0516
          TQMC-n2 : 0.0671
         ASTEROID : 0.9989
Model : 200 ntax, low ILS, deep speciation, 200 genes
            CA-ML : 0.0223
        ASTRID-ws : 0.0513
         ASTER-wh : 0.0352
       TQMC-wh_n2 : 0.0320
          TQMC-n2 : 0.0506
         ASTEROID : 0.9989
Model : 200 ntax, low ILS, deep speciation, 1000 genes
            CA-ML : 0.0178
        ASTRID-ws : 0.0485
         ASTER-wh : 0.0301
       TQMC-wh_n2 : 0.0248
          TQMC-n2 : 0.0452
         ASTEROID : 0.9989

Model : 200 ntax, low ILS, shallow speciation, 50 genes
            CA-ML : 0.0536
        ASTRID-ws : 0.0443
         ASTER-wh : 0.0481
       TQMC-wh_n2 : 0.0424
          TQMC-n2 : 0.0528
         ASTEROID : 0.9993
Model : 200 ntax, low ILS, shallow speciation, 200 genes
            CA-ML : 0.0311
        ASTRID-ws : 0.0223
         ASTER-wh : 0.0235
       TQMC-wh_n2 : 0.0224
          TQMC-n2 : 0.0291
         ASTEROID : 0.9995
Model : 200 ntax, low ILS, shallow speciation, 1000 genes
            CA-ML : 0.0144
        ASTRID-ws : 0.0130
         ASTER-wh : 0.0126
       TQMC-wh_n2 : 0.0112
          TQMC-n2 : 0.0188
         ASTEROID : 0.9995


Model : 200 ntax, medium ILS, deep speciation, 50 genes
            CA-ML : 0.1030
        ASTRID-ws : 0.0910
         ASTER-wh : 0.0892
       TQMC-wh_n2 : 0.0857
          TQMC-n2 : 0.0975
         ASTEROID : 0.9985
Model : 200 ntax, medium ILS, deep speciation, 200 genes
            CA-ML : 0.0569
        ASTRID-ws : 0.0550
         ASTER-wh : 0.0523
       TQMC-wh_n2 : 0.0466
          TQMC-n2 : 0.0570
         ASTEROID : 0.9987
Model : 200 ntax, medium ILS, deep speciation, 1000 genes
            CA-ML : 0.0284
        ASTRID-ws : 0.0392
         ASTER-wh : 0.0352
       TQMC-wh_n2 : 0.0309
          TQMC-n2 : 0.0398
         ASTEROID : 0.9987

Model : 200 ntax, medium ILS, shallow speciation, 50 genes
            CA-ML : 0.0922
        ASTRID-ws : 0.0728
         ASTER-wh : 0.0706
       TQMC-wh_n2 : 0.0669
          TQMC-n2 : 0.0809
         ASTEROID : 0.9985
Model : 200 ntax, medium ILS, shallow speciation, 200 genes
            CA-ML : 0.0552
        ASTRID-ws : 0.0417
         ASTER-wh : 0.0403
       TQMC-wh_n2 : 0.0382
          TQMC-n2 : 0.0475
         ASTEROID : 0.9987
Model : 200 ntax, medium ILS, shallow speciation, 1000 genes
            CA-ML : 0.0278
        ASTRID-ws : 0.0260
         ASTER-wh : 0.0241
       TQMC-wh_n2 : 0.0227
          TQMC-n2 : 0.0317
         ASTEROID : 0.9988


Model : 200 ntax, high ILS, deep speciation, 50 genes
            CA-ML : 0.2816
        ASTRID-ws : 0.2063
         ASTER-wh : 0.1898
       TQMC-wh_n2 : 0.1858
          TQMC-n2 : 0.2067
         ASTEROID : 0.9985
Model : 200 ntax, high ILS, deep speciation, 200 genes
            CA-ML : 0.1610
        ASTRID-ws : 0.1056
         ASTER-wh : 0.0913
       TQMC-wh_n2 : 0.0919
          TQMC-n2 : 0.1094
         ASTEROID : 0.9984
Model : 200 ntax, high ILS, deep speciation, 1000 genes
            CA-ML : 0.0801
        ASTRID-ws : 0.0493
         ASTER-wh : 0.0471
       TQMC-wh_n2 : 0.0486
          TQMC-n2 : 0.0617
         ASTEROID : 0.9986

Model : 200 ntax, high ILS, shallow speciation, 50 genes
            CA-ML : 0.2795
        ASTRID-ws : 0.2022
         ASTER-wh : 0.1818
       TQMC-wh_n2 : 0.1743
          TQMC-n2 : 0.1898
         ASTEROID : 0.9991
Model : 200 ntax, high ILS, shallow speciation, 200 genes
            CA-ML : 0.1625
        ASTRID-ws : 0.1013
         ASTER-wh : 0.0894
       TQMC-wh_n2 : 0.0830
          TQMC-n2 : 0.0960
         ASTEROID : 0.9989
Model : 200 ntax, high ILS, shallow speciation, 1000 genes
            CA-ML : 0.0795
        ASTRID-ws : 0.0471
         ASTER-wh : 0.0400
       TQMC-wh_n2 : 0.0377
          TQMC-n2 : 0.0493
         ASTEROID : 0.9990
"""