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
         "TQMC-ws_n2",
         "TQMC-n2",
         "TQMC-n2-origstudy"]

names = ["            CA-ML",
         "        ASTRID-ws",
         "         ASTER-wh",
         "       TQMC-wh_n2",
         "       TQMC-ws_n2",
         "          TQMC-n2",
         "TQMC-n2-origstudy"]

ngens = [50, 200, 1000]

for do in ["varyntax", "varyils"]: 
    if do == "varyntax":
        sys.stdout.write("Increasing number of taxa\n")
        ntaxs = [10, 50, 100, 200, 500, 1000]
        hghts = [2000000]
        rates = [0.000001]
        hghts = ["medium"]
        rates = ["shallow"]
    elif do == "varyils":
        sys.stdout.write("Increasing ILS\n")
        ntaxs = [200]
        hghts = ["low", "medium", "high"]
        rates = ["deep", "shallow"]

    df = pandas.read_csv("data-" + do + "-error.csv", keep_default_na=True)

    for ntax in ntaxs:
        for hght in hghts:
            for rate in rates:
                for ngen in ngens:
                    sys.stdout.write("Model : %s ntax, %s ILS, %s speciation, %d genes\n" % (ntax, hght, rate, ngen))
                    xdf = df[(df["NTAX"] == ntax) &
                             (df["ILSL"] == hght) &
                             (df["SPEC"] == rate) &
                             (df["NGEN"] == ngen)]

                    for mthd, name in zip(mthds, names):
                        if (mthd == "CA-ML") or \
                           (mthd == "TQMC-n2") or \
                           (mthd == "TQMC-n2-origstudy"):
                            supps = ["none"]
                            shift = "            "
                        else:
                            supps = ["sh", "abayes"]
                            shift = ""

                        sys.stdout.write("  %s : " % name)

                        nrepl = len(xdf[(xdf["MTHD"] == "CA-ML") &
                                        (xdf["SUPP"] == "none")].REPL.values)

                        for supp in supps:
                            if supp == "none":
                                ydf = xdf[(xdf["MTHD"] == mthd)]
                            else:
                                ydf = xdf[(xdf["MTHD"] == mthd) & (xdf["SUPP"] == supp)]

                            data = ydf.SEFNR.values
                            if nrepl != len(data) or nrepl > 50:
                                   sys.exit("wrong number of replicates")

                            with warnings.catch_warnings():
                                warnings.simplefilter("ignore", category=RuntimeWarning)
                                fnravg = numpy.nanmean(data)
                            fnravg = numpy.round(numpy.round(fnravg, 5), 4)

                            sys.stdout.write("%s%1.4f (%s) " % (shift, fnravg, supp))

                        sys.stdout.write("\n")
                    sys.stdout.write("\n")

"""
Increasing number of taxa
Model : 10 ntax, medium ILS, shallow speciation, 50 genes
              CA-ML :             0.0383 (none) 
          ASTRID-ws : 0.0306 (sh) 0.0281 (abayes) 
           ASTER-wh : 0.0255 (sh) 0.0306 (abayes) 
         TQMC-wh_n2 : 0.0230 (sh) 0.0306 (abayes) 
         TQMC-ws_n2 : 0.0230 (sh) 0.0281 (abayes) 
            TQMC-n2 :             0.0306 (none) 
  TQMC-n2-origstudy :             nan (none) 

Model : 10 ntax, medium ILS, shallow speciation, 200 genes
              CA-ML :             0.0179 (none) 
          ASTRID-ws : 0.0128 (sh) 0.0076 (abayes) 
           ASTER-wh : 0.0102 (sh) 0.0076 (abayes) 
         TQMC-wh_n2 : 0.0102 (sh) 0.0102 (abayes) 
         TQMC-ws_n2 : 0.0153 (sh) 0.0076 (abayes) 
            TQMC-n2 :             0.0153 (none) 
  TQMC-n2-origstudy :             nan (none) 

Model : 10 ntax, medium ILS, shallow speciation, 1000 genes
              CA-ML :             0.0208 (none) 
          ASTRID-ws : 0.0182 (sh) 0.0156 (abayes) 
           ASTER-wh : 0.0182 (sh) 0.0156 (abayes) 
         TQMC-wh_n2 : 0.0182 (sh) 0.0156 (abayes) 
         TQMC-ws_n2 : 0.0182 (sh) 0.0156 (abayes) 
            TQMC-n2 :             0.0156 (none) 
  TQMC-n2-origstudy :             0.0156 (none) 

Model : 50 ntax, medium ILS, shallow speciation, 50 genes
              CA-ML :             0.0780 (none) 
          ASTRID-ws : 0.0683 (sh) 0.0576 (abayes) 
           ASTER-wh : 0.0692 (sh) 0.0607 (abayes) 
         TQMC-wh_n2 : 0.0669 (sh) 0.0585 (abayes) 
         TQMC-ws_n2 : 0.0643 (sh) 0.0576 (abayes) 
            TQMC-n2 :             0.0705 (none) 
  TQMC-n2-origstudy :             nan (none) 

Model : 50 ntax, medium ILS, shallow speciation, 200 genes
              CA-ML :             0.0452 (none) 
          ASTRID-ws : 0.0394 (sh) 0.0355 (abayes) 
           ASTER-wh : 0.0399 (sh) 0.0306 (abayes) 
         TQMC-wh_n2 : 0.0399 (sh) 0.0284 (abayes) 
         TQMC-ws_n2 : 0.0390 (sh) 0.0315 (abayes) 
            TQMC-n2 :             0.0412 (none) 
  TQMC-n2-origstudy :             nan (none) 

Model : 50 ntax, medium ILS, shallow speciation, 1000 genes
              CA-ML :             0.0266 (none) 
          ASTRID-ws : 0.0230 (sh) 0.0191 (abayes) 
           ASTER-wh : 0.0235 (sh) 0.0186 (abayes) 
         TQMC-wh_n2 : 0.0226 (sh) 0.0177 (abayes) 
         TQMC-ws_n2 : 0.0208 (sh) 0.0177 (abayes) 
            TQMC-n2 :             0.0253 (none) 
  TQMC-n2-origstudy :             0.0248 (none) 

Model : 100 ntax, medium ILS, shallow speciation, 50 genes
              CA-ML :             0.0912 (none) 
          ASTRID-ws : 0.0742 (sh) 0.0704 (abayes) 
           ASTER-wh : 0.0757 (sh) 0.0716 (abayes) 
         TQMC-wh_n2 : 0.0678 (sh) 0.0672 (abayes) 
         TQMC-ws_n2 : 0.0663 (sh) 0.0638 (abayes) 
            TQMC-n2 :             0.0731 (none) 
  TQMC-n2-origstudy :             nan (none) 

Model : 100 ntax, medium ILS, shallow speciation, 200 genes
              CA-ML :             0.0474 (none) 
          ASTRID-ws : 0.0446 (sh) 0.0393 (abayes) 
           ASTER-wh : 0.0446 (sh) 0.0383 (abayes) 
         TQMC-wh_n2 : 0.0412 (sh) 0.0361 (abayes) 
         TQMC-ws_n2 : 0.0421 (sh) 0.0357 (abayes) 
            TQMC-n2 :             0.0457 (none) 
  TQMC-n2-origstudy :             nan (none) 

Model : 100 ntax, medium ILS, shallow speciation, 1000 genes
              CA-ML :             0.0249 (none) 
          ASTRID-ws : 0.0287 (sh) 0.0240 (abayes) 
           ASTER-wh : 0.0266 (sh) 0.0191 (abayes) 
         TQMC-wh_n2 : 0.0255 (sh) 0.0198 (abayes) 
         TQMC-ws_n2 : 0.0276 (sh) 0.0221 (abayes) 
            TQMC-n2 :             0.0304 (none) 
  TQMC-n2-origstudy :             0.0306 (none) 

Model : 200 ntax, medium ILS, shallow speciation, 50 genes
              CA-ML :             0.0922 (none) 
          ASTRID-ws : 0.0826 (sh) 0.0728 (abayes) 
           ASTER-wh : 0.0808 (sh) 0.0706 (abayes) 
         TQMC-wh_n2 : 0.0762 (sh) 0.0669 (abayes) 
         TQMC-ws_n2 : 0.0741 (sh) 0.0656 (abayes) 
            TQMC-n2 :             0.0802 (none) 
  TQMC-n2-origstudy :             nan (none) 

Model : 200 ntax, medium ILS, shallow speciation, 200 genes
              CA-ML :             0.0552 (none) 
          ASTRID-ws : 0.0465 (sh) 0.0417 (abayes) 
           ASTER-wh : 0.0473 (sh) 0.0403 (abayes) 
         TQMC-wh_n2 : 0.0431 (sh) 0.0382 (abayes) 
         TQMC-ws_n2 : 0.0427 (sh) 0.0369 (abayes) 
            TQMC-n2 :             0.0469 (none) 
  TQMC-n2-origstudy :             nan (none) 

Model : 200 ntax, medium ILS, shallow speciation, 1000 genes
              CA-ML :             0.0278 (none) 
          ASTRID-ws : 0.0304 (sh) 0.0260 (abayes) 
           ASTER-wh : 0.0303 (sh) 0.0241 (abayes) 
         TQMC-wh_n2 : 0.0262 (sh) 0.0227 (abayes) 
         TQMC-ws_n2 : 0.0279 (sh) 0.0236 (abayes) 
            TQMC-n2 :             0.0310 (none) 
  TQMC-n2-origstudy :             0.0315 (none) 

Model : 500 ntax, medium ILS, shallow speciation, 50 genes
              CA-ML :             0.0924 (none) 
          ASTRID-ws : 0.0821 (sh) 0.0740 (abayes) 
           ASTER-wh : 0.0789 (sh) 0.0740 (abayes) 
         TQMC-wh_n2 : 0.0723 (sh) 0.0665 (abayes) 
         TQMC-ws_n2 : 0.0712 (sh) 0.0664 (abayes) 
            TQMC-n2 :             0.0752 (none) 
  TQMC-n2-origstudy :             nan (none) 

Model : 500 ntax, medium ILS, shallow speciation, 200 genes
              CA-ML :             0.0472 (none) 
          ASTRID-ws : 0.0457 (sh) 0.0410 (abayes) 
           ASTER-wh : 0.0442 (sh) 0.0400 (abayes) 
         TQMC-wh_n2 : 0.0378 (sh) 0.0341 (abayes) 
         TQMC-ws_n2 : 0.0393 (sh) 0.0374 (abayes) 
            TQMC-n2 :             0.0422 (none) 
  TQMC-n2-origstudy :             nan (none) 

Model : 500 ntax, medium ILS, shallow speciation, 1000 genes
              CA-ML :             0.0233 (none) 
          ASTRID-ws : 0.0300 (sh) 0.0261 (abayes) 
           ASTER-wh : 0.0276 (sh) 0.0243 (abayes) 
         TQMC-wh_n2 : 0.0227 (sh) 0.0203 (abayes) 
         TQMC-ws_n2 : 0.0258 (sh) 0.0237 (abayes) 
            TQMC-n2 :             0.0276 (none) 
  TQMC-n2-origstudy :             0.0280 (none) 

Model : 1000 ntax, medium ILS, shallow speciation, 50 genes
              CA-ML :             0.0976 (none) 
          ASTRID-ws : 0.1001 (sh) 0.0881 (abayes) 
           ASTER-wh : 0.0979 (sh) 0.0867 (abayes) 
         TQMC-wh_n2 : 0.0884 (sh) 0.0768 (abayes) 
         TQMC-ws_n2 : 0.0912 (sh) 0.0804 (abayes) 
            TQMC-n2 :             0.0985 (none) 
  TQMC-n2-origstudy :             nan (none) 

Model : 1000 ntax, medium ILS, shallow speciation, 200 genes
              CA-ML :             0.0515 (none) 
          ASTRID-ws : 0.0588 (sh) 0.0518 (abayes) 
           ASTER-wh : 0.0557 (sh) 0.0485 (abayes) 
         TQMC-wh_n2 : 0.0488 (sh) 0.0418 (abayes) 
         TQMC-ws_n2 : 0.0537 (sh) 0.0469 (abayes) 
            TQMC-n2 :             0.0576 (none) 
  TQMC-n2-origstudy :             nan (none) 

Model : 1000 ntax, medium ILS, shallow speciation, 1000 genes
              CA-ML :             nan (none) 
          ASTRID-ws : 0.0403 (sh) 0.0351 (abayes) 
           ASTER-wh : 0.0375 (sh) 0.0321 (abayes) 
         TQMC-wh_n2 : 0.0315 (sh) 0.0258 (abayes) 
         TQMC-ws_n2 : 0.0369 (sh) 0.0318 (abayes) 
            TQMC-n2 :             0.0411 (none) 
  TQMC-n2-origstudy :             0.0416 (none) 

Increasing ILS
Model : 200 ntax, low ILS, deep speciation, 50 genes
              CA-ML :             0.0398 (none) 
          ASTRID-ws : 0.0748 (sh) 0.0649 (abayes) 
           ASTER-wh : 0.0676 (sh) 0.0583 (abayes) 
         TQMC-wh_n2 : 0.0628 (sh) 0.0516 (abayes) 
         TQMC-ws_n2 : 0.0654 (sh) 0.0586 (abayes) 
            TQMC-n2 :             0.0656 (none) 
  TQMC-n2-origstudy :             nan (none) 

Model : 200 ntax, low ILS, deep speciation, 200 genes
              CA-ML :             0.0223 (none) 
          ASTRID-ws : 0.0568 (sh) 0.0513 (abayes) 
           ASTER-wh : 0.0435 (sh) 0.0352 (abayes) 
         TQMC-wh_n2 : 0.0393 (sh) 0.0320 (abayes) 
         TQMC-ws_n2 : 0.0471 (sh) 0.0432 (abayes) 
            TQMC-n2 :             0.0505 (none) 
  TQMC-n2-origstudy :             nan (none) 

Model : 200 ntax, low ILS, deep speciation, 1000 genes
              CA-ML :             0.0178 (none) 
          ASTRID-ws : 0.0538 (sh) 0.0485 (abayes) 
           ASTER-wh : 0.0363 (sh) 0.0301 (abayes) 
         TQMC-wh_n2 : 0.0332 (sh) 0.0248 (abayes) 
         TQMC-ws_n2 : 0.0439 (sh) 0.0401 (abayes) 
            TQMC-n2 :             0.0446 (none) 
  TQMC-n2-origstudy :             0.0450 (none) 

Model : 200 ntax, low ILS, shallow speciation, 50 genes
              CA-ML :             0.0536 (none) 
          ASTRID-ws : 0.0507 (sh) 0.0443 (abayes) 
           ASTER-wh : 0.0526 (sh) 0.0481 (abayes) 
         TQMC-wh_n2 : 0.0481 (sh) 0.0424 (abayes) 
         TQMC-ws_n2 : 0.0470 (sh) 0.0420 (abayes) 
            TQMC-n2 :             0.0502 (none) 
  TQMC-n2-origstudy :             nan (none) 

Model : 200 ntax, low ILS, shallow speciation, 200 genes
              CA-ML :             0.0311 (none) 
          ASTRID-ws : 0.0273 (sh) 0.0223 (abayes) 
           ASTER-wh : 0.0295 (sh) 0.0235 (abayes) 
         TQMC-wh_n2 : 0.0265 (sh) 0.0224 (abayes) 
         TQMC-ws_n2 : 0.0261 (sh) 0.0213 (abayes) 
            TQMC-n2 :             0.0276 (none) 
  TQMC-n2-origstudy :             nan (none) 

Model : 200 ntax, low ILS, shallow speciation, 1000 genes
              CA-ML :             0.0144 (none) 
          ASTRID-ws : 0.0160 (sh) 0.0130 (abayes) 
           ASTER-wh : 0.0176 (sh) 0.0126 (abayes) 
         TQMC-wh_n2 : 0.0148 (sh) 0.0112 (abayes) 
         TQMC-ws_n2 : 0.0157 (sh) 0.0120 (abayes) 
            TQMC-n2 :             0.0175 (none) 
  TQMC-n2-origstudy :             0.0183 (none) 

Model : 200 ntax, medium ILS, deep speciation, 50 genes
              CA-ML :             0.1030 (none) 
          ASTRID-ws : 0.1014 (sh) 0.0910 (abayes) 
           ASTER-wh : 0.0985 (sh) 0.0892 (abayes) 
         TQMC-wh_n2 : 0.0925 (sh) 0.0857 (abayes) 
         TQMC-ws_n2 : 0.0911 (sh) 0.0834 (abayes) 
            TQMC-n2 :             0.0962 (none) 
  TQMC-n2-origstudy :             nan (none) 

Model : 200 ntax, medium ILS, deep speciation, 200 genes
              CA-ML :             0.0569 (none) 
          ASTRID-ws : 0.0609 (sh) 0.0550 (abayes) 
           ASTER-wh : 0.0561 (sh) 0.0523 (abayes) 
         TQMC-wh_n2 : 0.0506 (sh) 0.0466 (abayes) 
         TQMC-ws_n2 : 0.0527 (sh) 0.0488 (abayes) 
            TQMC-n2 :             0.0571 (none) 
  TQMC-n2-origstudy :             nan (none) 

Model : 200 ntax, medium ILS, deep speciation, 1000 genes
              CA-ML :             0.0284 (none) 
          ASTRID-ws : 0.0434 (sh) 0.0392 (abayes) 
           ASTER-wh : 0.0377 (sh) 0.0352 (abayes) 
         TQMC-wh_n2 : 0.0341 (sh) 0.0309 (abayes) 
         TQMC-ws_n2 : 0.0377 (sh) 0.0357 (abayes) 
            TQMC-n2 :             0.0397 (none) 
  TQMC-n2-origstudy :             0.0396 (none) 

Model : 200 ntax, medium ILS, shallow speciation, 50 genes
              CA-ML :             0.0922 (none) 
          ASTRID-ws : 0.0826 (sh) 0.0728 (abayes) 
           ASTER-wh : 0.0808 (sh) 0.0706 (abayes) 
         TQMC-wh_n2 : 0.0762 (sh) 0.0669 (abayes) 
         TQMC-ws_n2 : 0.0741 (sh) 0.0656 (abayes) 
            TQMC-n2 :             0.0802 (none) 
  TQMC-n2-origstudy :             nan (none) 

Model : 200 ntax, medium ILS, shallow speciation, 200 genes
              CA-ML :             0.0552 (none) 
          ASTRID-ws : 0.0465 (sh) 0.0417 (abayes) 
           ASTER-wh : 0.0473 (sh) 0.0403 (abayes) 
         TQMC-wh_n2 : 0.0431 (sh) 0.0382 (abayes) 
         TQMC-ws_n2 : 0.0427 (sh) 0.0369 (abayes) 
            TQMC-n2 :             0.0469 (none) 
  TQMC-n2-origstudy :             nan (none) 

Model : 200 ntax, medium ILS, shallow speciation, 1000 genes
              CA-ML :             0.0278 (none) 
          ASTRID-ws : 0.0304 (sh) 0.0260 (abayes) 
           ASTER-wh : 0.0303 (sh) 0.0241 (abayes) 
         TQMC-wh_n2 : 0.0262 (sh) 0.0227 (abayes) 
         TQMC-ws_n2 : 0.0279 (sh) 0.0236 (abayes) 
            TQMC-n2 :             0.0310 (none) 
  TQMC-n2-origstudy :             0.0315 (none) 

Model : 200 ntax, high ILS, deep speciation, 50 genes
              CA-ML :             0.2816 (none) 
          ASTRID-ws : 0.2154 (sh) 0.2063 (abayes) 
           ASTER-wh : 0.2019 (sh) 0.1898 (abayes) 
         TQMC-wh_n2 : 0.1956 (sh) 0.1858 (abayes) 
         TQMC-ws_n2 : 0.1975 (sh) 0.1833 (abayes) 
            TQMC-n2 :             0.2042 (none) 
  TQMC-n2-origstudy :             nan (none) 

Model : 200 ntax, high ILS, deep speciation, 200 genes
              CA-ML :             0.1610 (none) 
          ASTRID-ws : 0.1133 (sh) 0.1056 (abayes) 
           ASTER-wh : 0.1016 (sh) 0.0913 (abayes) 
         TQMC-wh_n2 : 0.0981 (sh) 0.0919 (abayes) 
         TQMC-ws_n2 : 0.0988 (sh) 0.0915 (abayes) 
            TQMC-n2 :             0.1076 (none) 
  TQMC-n2-origstudy :             nan (none) 

Model : 200 ntax, high ILS, deep speciation, 1000 genes
              CA-ML :             0.0799 (none) 
          ASTRID-ws : 0.0567 (sh) 0.0494 (abayes) 
           ASTER-wh : 0.0554 (sh) 0.0473 (abayes) 
         TQMC-wh_n2 : 0.0554 (sh) 0.0487 (abayes) 
         TQMC-ws_n2 : 0.0562 (sh) 0.0492 (abayes) 
            TQMC-n2 :             0.0609 (none) 
  TQMC-n2-origstudy :             0.0628 (none) 

Model : 200 ntax, high ILS, shallow speciation, 50 genes
              CA-ML :             0.2795 (none) 
          ASTRID-ws : 0.2092 (sh) 0.2022 (abayes) 
           ASTER-wh : 0.1880 (sh) 0.1818 (abayes) 
         TQMC-wh_n2 : 0.1785 (sh) 0.1743 (abayes) 
         TQMC-ws_n2 : 0.1789 (sh) 0.1706 (abayes) 
            TQMC-n2 :             0.1839 (none) 
  TQMC-n2-origstudy :             nan (none) 

Model : 200 ntax, high ILS, shallow speciation, 200 genes
              CA-ML :             0.1625 (none) 
          ASTRID-ws : 0.1094 (sh) 0.1013 (abayes) 
           ASTER-wh : 0.0960 (sh) 0.0894 (abayes) 
         TQMC-wh_n2 : 0.0916 (sh) 0.0830 (abayes) 
         TQMC-ws_n2 : 0.0888 (sh) 0.0811 (abayes) 
            TQMC-n2 :             0.0965 (none) 
  TQMC-n2-origstudy :             nan (none) 

Model : 200 ntax, high ILS, shallow speciation, 1000 genes
              CA-ML :             0.0795 (none) 
          ASTRID-ws : 0.0512 (sh) 0.0471 (abayes) 
           ASTER-wh : 0.0471 (sh) 0.0400 (abayes) 
         TQMC-wh_n2 : 0.0436 (sh) 0.0377 (abayes) 
         TQMC-ws_n2 : 0.0459 (sh) 0.0404 (abayes) 
            TQMC-n2 :             0.0493 (none) 
  TQMC-n2-origstudy :             0.0498 (none)
"""