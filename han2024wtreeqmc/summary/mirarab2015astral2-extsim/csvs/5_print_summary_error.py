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

mthds = ["asteroid",
         "wastrid_s",
         "aster_h_t16",
         "wtreeqmc_wh_n2",
         "caml",
         "wtreeqmc_wh_n1",
         "wtreeqmc_wh_n0",
         "wtreeqmc_ws_n2",
         "wtreeqmc_wn_n2"]

names = ["       asteroid",
         "      wastrid_s",
         "   *aster_h_t16",
         "*wtreeqmc_wh_n2",
         "          *caml",
         " wtreeqmc_wh_n1",
         " wtreeqmc_wh_n0",
         " wtreeqmc_ws_n2",
         " wtreeqmc_wn_n2"]

ngens = [50, 200, 1000]

for do in ["ntax", "ils"]: 
    if do == "ntax":
        sys.stdout.write("Increasing number of taxa\n")
        df = pandas.read_csv("data-varyntax-error.csv", keep_default_na=True)
        ntaxs = [10, 50, 100, 200, 500, 1000]
        hghts = [2000000]
        rates = [0.000001]
        
    elif do == "ils":
        sys.stdout.write("Increasing ILS\n")
        df = pandas.read_csv("data-varyils-error.csv", keep_default_na=True)
        ntaxs = [200]
        hghts = [10000000, 2000000, 500000]
        rates = [0.0000001, 0.000001]

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
                        if (mthd == "asteroid") or (mthd == "wtreeqmc_wn_n2") or (mthd == "caml"):
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
Model : 10 ntax, 2000000 hght, 1e-06 rate, 50 genes
         asteroid :             0.4133 (none) 
        wastrid_s : 0.0306 (sh) 0.0281 (abayes) 
     *aster_h_t16 : 0.0255 (sh) 0.0306 (abayes) 
  *wtreeqmc_wh_n2 : 0.0230 (sh) 0.0306 (abayes) 
            *caml :             0.0383 (none) 
   wtreeqmc_wh_n1 : 0.0230 (sh) 0.0281 (abayes) 
   wtreeqmc_wh_n0 : 0.0230 (sh) 0.0281 (abayes) 
   wtreeqmc_ws_n2 : 0.0230 (sh) 0.0281 (abayes) 
   wtreeqmc_wn_n2 :             0.0306 (none) 

Model : 10 ntax, 2000000 hght, 1e-06 rate, 200 genes
         asteroid :             0.4082 (none) 
        wastrid_s : 0.0128 (sh) 0.0076 (abayes) 
     *aster_h_t16 : 0.0102 (sh) 0.0076 (abayes) 
  *wtreeqmc_wh_n2 : 0.0102 (sh) 0.0102 (abayes) 
            *caml :             0.0179 (none) 
   wtreeqmc_wh_n1 : 0.0153 (sh) 0.0076 (abayes) 
   wtreeqmc_wh_n0 : 0.0153 (sh) 0.0076 (abayes) 
   wtreeqmc_ws_n2 : 0.0153 (sh) 0.0076 (abayes) 
   wtreeqmc_wn_n2 :             0.0153 (none) 

Model : 10 ntax, 2000000 hght, 1e-06 rate, 1000 genes
         asteroid :             0.4115 (none) 
        wastrid_s : 0.0182 (sh) 0.0156 (abayes) 
     *aster_h_t16 : 0.0182 (sh) 0.0156 (abayes) 
  *wtreeqmc_wh_n2 : 0.0182 (sh) 0.0156 (abayes) 
            *caml :             0.0208 (none) 
   wtreeqmc_wh_n1 : 0.0182 (sh) 0.0156 (abayes) 
   wtreeqmc_wh_n0 : 0.0182 (sh) 0.0156 (abayes) 
   wtreeqmc_ws_n2 : 0.0182 (sh) 0.0156 (abayes) 
   wtreeqmc_wn_n2 :             0.0156 (none) 

Model : 50 ntax, 2000000 hght, 1e-06 rate, 50 genes
         asteroid :             0.9942 (none) 
        wastrid_s : 0.0683 (sh) 0.0576 (abayes) 
     *aster_h_t16 : 0.0692 (sh) 0.0607 (abayes) 
  *wtreeqmc_wh_n2 : 0.0669 (sh) 0.0585 (abayes) 
            *caml :             0.0780 (none) 
   wtreeqmc_wh_n1 : 0.0643 (sh) 0.0576 (abayes) 
   wtreeqmc_wh_n0 : 0.0643 (sh) 0.0576 (abayes) 
   wtreeqmc_ws_n2 : 0.0643 (sh) 0.0576 (abayes) 
   wtreeqmc_wn_n2 :             0.0705 (none) 

Model : 50 ntax, 2000000 hght, 1e-06 rate, 200 genes
         asteroid :             0.9934 (none) 
        wastrid_s : 0.0394 (sh) 0.0355 (abayes) 
     *aster_h_t16 : 0.0399 (sh) 0.0306 (abayes) 
  *wtreeqmc_wh_n2 : 0.0399 (sh) 0.0284 (abayes) 
            *caml :             0.0452 (none) 
   wtreeqmc_wh_n1 : 0.0390 (sh) 0.0315 (abayes) 
   wtreeqmc_wh_n0 : 0.0390 (sh) 0.0315 (abayes) 
   wtreeqmc_ws_n2 : 0.0390 (sh) 0.0315 (abayes) 
   wtreeqmc_wn_n2 :             0.0412 (none) 

Model : 50 ntax, 2000000 hght, 1e-06 rate, 1000 genes
         asteroid :             0.9934 (none) 
        wastrid_s : 0.0230 (sh) 0.0191 (abayes) 
     *aster_h_t16 : 0.0235 (sh) 0.0186 (abayes) 
  *wtreeqmc_wh_n2 : 0.0226 (sh) 0.0177 (abayes) 
            *caml :             0.0266 (none) 
   wtreeqmc_wh_n1 : 0.0208 (sh) 0.0177 (abayes) 
   wtreeqmc_wh_n0 : 0.0208 (sh) 0.0177 (abayes) 
   wtreeqmc_ws_n2 : 0.0208 (sh) 0.0177 (abayes) 
   wtreeqmc_wn_n2 :             0.0253 (none) 

Model : 100 ntax, 2000000 hght, 1e-06 rate, 50 genes
         asteroid :             0.9964 (none) 
        wastrid_s : 0.0742 (sh) 0.0704 (abayes) 
     *aster_h_t16 : 0.0757 (sh) 0.0716 (abayes) 
  *wtreeqmc_wh_n2 : 0.0678 (sh) 0.0672 (abayes) 
            *caml :             0.0912 (none) 
   wtreeqmc_wh_n1 : 0.0663 (sh) 0.0638 (abayes) 
   wtreeqmc_wh_n0 : 0.0663 (sh) 0.0638 (abayes) 
   wtreeqmc_ws_n2 : 0.0663 (sh) 0.0638 (abayes) 
   wtreeqmc_wn_n2 :             0.0731 (none) 

Model : 100 ntax, 2000000 hght, 1e-06 rate, 200 genes
         asteroid :             0.9966 (none) 
        wastrid_s : 0.0446 (sh) 0.0393 (abayes) 
     *aster_h_t16 : 0.0446 (sh) 0.0383 (abayes) 
  *wtreeqmc_wh_n2 : 0.0412 (sh) 0.0361 (abayes) 
            *caml :             0.0474 (none) 
   wtreeqmc_wh_n1 : 0.0421 (sh) 0.0357 (abayes) 
   wtreeqmc_wh_n0 : 0.0421 (sh) 0.0357 (abayes) 
   wtreeqmc_ws_n2 : 0.0421 (sh) 0.0357 (abayes) 
   wtreeqmc_wn_n2 :             0.0457 (none) 

Model : 100 ntax, 2000000 hght, 1e-06 rate, 1000 genes
         asteroid :             0.9964 (none) 
        wastrid_s : 0.0287 (sh) 0.0240 (abayes) 
     *aster_h_t16 : 0.0266 (sh) 0.0191 (abayes) 
  *wtreeqmc_wh_n2 : 0.0255 (sh) 0.0198 (abayes) 
            *caml :             0.0249 (none) 
   wtreeqmc_wh_n1 : 0.0276 (sh) 0.0221 (abayes) 
   wtreeqmc_wh_n0 : 0.0276 (sh) 0.0221 (abayes) 
   wtreeqmc_ws_n2 : 0.0276 (sh) 0.0221 (abayes) 
   wtreeqmc_wn_n2 :             0.0304 (none) 

Model : 200 ntax, 2000000 hght, 1e-06 rate, 50 genes
         asteroid :             0.9984 (none) 
        wastrid_s : 0.0826 (sh) 0.0728 (abayes) 
     *aster_h_t16 : 0.0808 (sh) 0.0706 (abayes) 
  *wtreeqmc_wh_n2 : 0.0762 (sh) 0.0669 (abayes) 
            *caml :             0.0922 (none) 
   wtreeqmc_wh_n1 : 0.0741 (sh) 0.0656 (abayes) 
   wtreeqmc_wh_n0 : 0.0741 (sh) 0.0656 (abayes) 
   wtreeqmc_ws_n2 : 0.0741 (sh) 0.0656 (abayes) 
   wtreeqmc_wn_n2 :             0.0802 (none) 

Model : 200 ntax, 2000000 hght, 1e-06 rate, 200 genes
         asteroid :             0.9986 (none) 
        wastrid_s : 0.0465 (sh) 0.0417 (abayes) 
     *aster_h_t16 : 0.0473 (sh) 0.0403 (abayes) 
  *wtreeqmc_wh_n2 : 0.0431 (sh) 0.0382 (abayes) 
            *caml :             0.0552 (none) 
   wtreeqmc_wh_n1 : 0.0427 (sh) 0.0369 (abayes) 
   wtreeqmc_wh_n0 : 0.0427 (sh) 0.0369 (abayes) 
   wtreeqmc_ws_n2 : 0.0427 (sh) 0.0369 (abayes) 
   wtreeqmc_wn_n2 :             0.0469 (none) 

Model : 200 ntax, 2000000 hght, 1e-06 rate, 1000 genes
         asteroid :             0.9987 (none) 
        wastrid_s : 0.0304 (sh) 0.0260 (abayes) 
     *aster_h_t16 : 0.0303 (sh) 0.0241 (abayes) 
  *wtreeqmc_wh_n2 : 0.0262 (sh) 0.0227 (abayes) 
            *caml :             0.0278 (none) 
   wtreeqmc_wh_n1 : 0.0279 (sh) 0.0236 (abayes) 
   wtreeqmc_wh_n0 : 0.0279 (sh) 0.0236 (abayes) 
   wtreeqmc_ws_n2 : 0.0279 (sh) 0.0236 (abayes) 
   wtreeqmc_wn_n2 :             0.0310 (none) 

Model : 500 ntax, 2000000 hght, 1e-06 rate, 50 genes
         asteroid :             0.9995 (none) 
        wastrid_s : 0.0821 (sh) 0.0740 (abayes) 
     *aster_h_t16 : 0.0789 (sh) 0.0740 (abayes) 
  *wtreeqmc_wh_n2 : 0.0723 (sh) 0.0665 (abayes) 
            *caml :             0.0924 (none) 
   wtreeqmc_wh_n1 : 0.0712 (sh) 0.0664 (abayes) 
   wtreeqmc_wh_n0 : 0.0712 (sh) 0.0664 (abayes) 
   wtreeqmc_ws_n2 : 0.0712 (sh) 0.0664 (abayes) 
   wtreeqmc_wn_n2 :             0.0752 (none) 

Model : 500 ntax, 2000000 hght, 1e-06 rate, 200 genes
         asteroid :             0.9996 (none) 
        wastrid_s : 0.0457 (sh) 0.0410 (abayes) 
     *aster_h_t16 : 0.0442 (sh) 0.0400 (abayes) 
  *wtreeqmc_wh_n2 : 0.0378 (sh) 0.0341 (abayes) 
            *caml :             0.0472 (none) 
   wtreeqmc_wh_n1 : 0.0393 (sh) 0.0374 (abayes) 
   wtreeqmc_wh_n0 : 0.0393 (sh) 0.0374 (abayes) 
   wtreeqmc_ws_n2 : 0.0393 (sh) 0.0374 (abayes) 
   wtreeqmc_wn_n2 :             0.0422 (none) 

Model : 500 ntax, 2000000 hght, 1e-06 rate, 1000 genes
         asteroid :             0.9996 (none) 
        wastrid_s : 0.0300 (sh) 0.0261 (abayes) 
     *aster_h_t16 : 0.0276 (sh) 0.0243 (abayes) 
  *wtreeqmc_wh_n2 : 0.0227 (sh) 0.0203 (abayes) 
            *caml :             0.0233 (none) 
   wtreeqmc_wh_n1 : 0.0258 (sh) 0.0237 (abayes) 
   wtreeqmc_wh_n0 : 0.0258 (sh) 0.0237 (abayes) 
   wtreeqmc_ws_n2 : 0.0258 (sh) 0.0237 (abayes) 
   wtreeqmc_wn_n2 :             0.0276 (none) 

Model : 1000 ntax, 2000000 hght, 1e-06 rate, 50 genes
         asteroid :             0.9998 (none) 
        wastrid_s : 0.1001 (sh) 0.0881 (abayes) 
     *aster_h_t16 : 0.0979 (sh) 0.0867 (abayes) 
  *wtreeqmc_wh_n2 : 0.0884 (sh) 0.0768 (abayes) 
            *caml :             0.0976 (none) 
   wtreeqmc_wh_n1 : 0.0912 (sh) 0.0804 (abayes) 
   wtreeqmc_wh_n0 : 0.0912 (sh) 0.0804 (abayes) 
   wtreeqmc_ws_n2 : 0.0912 (sh) 0.0804 (abayes) 
   wtreeqmc_wn_n2 :             0.0985 (none) 

Model : 1000 ntax, 2000000 hght, 1e-06 rate, 200 genes
         asteroid :             0.9998 (none) 
        wastrid_s : 0.0588 (sh) 0.0518 (abayes) 
     *aster_h_t16 : 0.0557 (sh) 0.0485 (abayes) 
  *wtreeqmc_wh_n2 : 0.0488 (sh) 0.0418 (abayes) 
            *caml :             0.0515 (none) 
   wtreeqmc_wh_n1 : 0.0537 (sh) 0.0469 (abayes) 
   wtreeqmc_wh_n0 : 0.0537 (sh) 0.0469 (abayes) 
   wtreeqmc_ws_n2 : 0.0537 (sh) 0.0469 (abayes) 
   wtreeqmc_wn_n2 :             0.0576 (none) 

Model : 1000 ntax, 2000000 hght, 1e-06 rate, 1000 genes
         asteroid :             0.9998 (none) 
        wastrid_s : 0.0403 (sh) 0.0351 (abayes) 
     *aster_h_t16 : 0.0375 (sh) 0.0321 (abayes) 
  *wtreeqmc_wh_n2 : 0.0315 (sh) 0.0258 (abayes)
            *caml :             nan (none) 
   wtreeqmc_wh_n1 : 0.0369 (sh) 0.0318 (abayes) 
   wtreeqmc_wh_n0 : 0.0369 (sh) 0.0318 (abayes) 
   wtreeqmc_ws_n2 : 0.0369 (sh) 0.0318 (abayes) 
   wtreeqmc_wn_n2 :             0.0411 (none) 

Increasing ILS
Model : 200 ntax, 10000000 hght, 1e-07 rate, 50 genes
         asteroid :             0.9989 (none) 
        wastrid_s : 0.0748 (sh) 0.0649 (abayes) 
     *aster_h_t16 : 0.0676 (sh) 0.0583 (abayes) 
  *wtreeqmc_wh_n2 : 0.0628 (sh) 0.0516 (abayes) 
            *caml :             0.0398 (none) 
   wtreeqmc_wh_n1 : 0.0654 (sh) 0.0586 (abayes) 
   wtreeqmc_wh_n0 : 0.0654 (sh) 0.0586 (abayes) 
   wtreeqmc_ws_n2 : 0.0654 (sh) 0.0586 (abayes) 
   wtreeqmc_wn_n2 :             0.0656 (none) 

Model : 200 ntax, 10000000 hght, 1e-07 rate, 200 genes
         asteroid :             0.9989 (none) 
        wastrid_s : 0.0568 (sh) 0.0513 (abayes) 
     *aster_h_t16 : 0.0435 (sh) 0.0352 (abayes) 
  *wtreeqmc_wh_n2 : 0.0393 (sh) 0.0320 (abayes) 
            *caml :             0.0223 (none) 
   wtreeqmc_wh_n1 : 0.0471 (sh) 0.0432 (abayes) 
   wtreeqmc_wh_n0 : 0.0471 (sh) 0.0432 (abayes) 
   wtreeqmc_ws_n2 : 0.0471 (sh) 0.0432 (abayes) 
   wtreeqmc_wn_n2 :             0.0505 (none) 

Model : 200 ntax, 10000000 hght, 1e-07 rate, 1000 genes
         asteroid :             0.9989 (none) 
        wastrid_s : 0.0538 (sh) 0.0485 (abayes) 
     *aster_h_t16 : 0.0363 (sh) 0.0301 (abayes) 
  *wtreeqmc_wh_n2 : 0.0332 (sh) 0.0248 (abayes) 
            *caml :             0.0178 (none) 
   wtreeqmc_wh_n1 : 0.0439 (sh) 0.0401 (abayes) 
   wtreeqmc_wh_n0 : 0.0439 (sh) 0.0401 (abayes) 
   wtreeqmc_ws_n2 : 0.0439 (sh) 0.0401 (abayes) 
   wtreeqmc_wn_n2 :             0.0446 (none) 

Model : 200 ntax, 10000000 hght, 1e-06 rate, 50 genes
         asteroid :             0.9994 (none) 
        wastrid_s : 0.0507 (sh) 0.0443 (abayes) 
     *aster_h_t16 : 0.0526 (sh) 0.0481 (abayes) 
  *wtreeqmc_wh_n2 : 0.0481 (sh) 0.0424 (abayes) 
            *caml :             0.0536 (none) 
   wtreeqmc_wh_n1 : 0.0470 (sh) 0.0420 (abayes) 
   wtreeqmc_wh_n0 : 0.0470 (sh) 0.0420 (abayes) 
   wtreeqmc_ws_n2 : 0.0470 (sh) 0.0420 (abayes) 
   wtreeqmc_wn_n2 :             0.0502 (none) 

Model : 200 ntax, 10000000 hght, 1e-06 rate, 200 genes
         asteroid :             0.9995 (none) 
        wastrid_s : 0.0273 (sh) 0.0223 (abayes) 
     *aster_h_t16 : 0.0295 (sh) 0.0235 (abayes) 
  *wtreeqmc_wh_n2 : 0.0265 (sh) 0.0224 (abayes) 
            *caml :             0.0311 (none) 
   wtreeqmc_wh_n1 : 0.0261 (sh) 0.0213 (abayes) 
   wtreeqmc_wh_n0 : 0.0261 (sh) 0.0213 (abayes) 
   wtreeqmc_ws_n2 : 0.0261 (sh) 0.0213 (abayes) 
   wtreeqmc_wn_n2 :             0.0276 (none) 

Model : 200 ntax, 10000000 hght, 1e-06 rate, 1000 genes
         asteroid :             0.9995 (none) 
        wastrid_s : 0.0160 (sh) 0.0130 (abayes) 
     *aster_h_t16 : 0.0176 (sh) 0.0126 (abayes) 
  *wtreeqmc_wh_n2 : 0.0148 (sh) 0.0112 (abayes) 
            *caml :             0.0144 (none) 
   wtreeqmc_wh_n1 : 0.0157 (sh) 0.0120 (abayes) 
   wtreeqmc_wh_n0 : 0.0157 (sh) 0.0120 (abayes) 
   wtreeqmc_ws_n2 : 0.0157 (sh) 0.0120 (abayes) 
   wtreeqmc_wn_n2 :             0.0175 (none) 

Model : 200 ntax, 2000000 hght, 1e-07 rate, 50 genes
         asteroid :             0.9985 (none) 
        wastrid_s : 0.1014 (sh) 0.0910 (abayes) 
     *aster_h_t16 : 0.0985 (sh) 0.0892 (abayes) 
  *wtreeqmc_wh_n2 : 0.0925 (sh) 0.0857 (abayes) 
            *caml :             0.1030 (none) 
   wtreeqmc_wh_n1 : 0.0911 (sh) 0.0834 (abayes) 
   wtreeqmc_wh_n0 : 0.0911 (sh) 0.0834 (abayes) 
   wtreeqmc_ws_n2 : 0.0911 (sh) 0.0834 (abayes) 
   wtreeqmc_wn_n2 :             0.0962 (none) 

Model : 200 ntax, 2000000 hght, 1e-07 rate, 200 genes
         asteroid :             0.9985 (none) 
        wastrid_s : 0.0609 (sh) 0.0550 (abayes) 
     *aster_h_t16 : 0.0561 (sh) 0.0523 (abayes) 
  *wtreeqmc_wh_n2 : 0.0506 (sh) 0.0466 (abayes) 
            *caml :             0.0569 (none) 
   wtreeqmc_wh_n1 : 0.0527 (sh) 0.0488 (abayes) 
   wtreeqmc_wh_n0 : 0.0527 (sh) 0.0488 (abayes) 
   wtreeqmc_ws_n2 : 0.0527 (sh) 0.0488 (abayes) 
   wtreeqmc_wn_n2 :             0.0571 (none) 

Model : 200 ntax, 2000000 hght, 1e-07 rate, 1000 genes
         asteroid :             0.9986 (none) 
        wastrid_s : 0.0434 (sh) 0.0392 (abayes) 
     *aster_h_t16 : 0.0377 (sh) 0.0352 (abayes) 
  *wtreeqmc_wh_n2 : 0.0341 (sh) 0.0309 (abayes) 
            *caml :             0.0284 (none) 
   wtreeqmc_wh_n1 : 0.0377 (sh) 0.0357 (abayes) 
   wtreeqmc_wh_n0 : 0.0377 (sh) 0.0357 (abayes) 
   wtreeqmc_ws_n2 : 0.0377 (sh) 0.0357 (abayes) 
   wtreeqmc_wn_n2 :             0.0397 (none) 

Model : 200 ntax, 2000000 hght, 1e-06 rate, 50 genes
         asteroid :             0.9984 (none) 
        wastrid_s : 0.0826 (sh) 0.0728 (abayes) 
     *aster_h_t16 : 0.0808 (sh) 0.0706 (abayes) 
  *wtreeqmc_wh_n2 : 0.0762 (sh) 0.0669 (abayes) 
            *caml :             0.0922 (none) 
   wtreeqmc_wh_n1 : 0.0741 (sh) 0.0656 (abayes) 
   wtreeqmc_wh_n0 : 0.0741 (sh) 0.0656 (abayes) 
   wtreeqmc_ws_n2 : 0.0741 (sh) 0.0656 (abayes) 
   wtreeqmc_wn_n2 :             0.0802 (none) 

Model : 200 ntax, 2000000 hght, 1e-06 rate, 200 genes
         asteroid :             0.9986 (none) 
        wastrid_s : 0.0465 (sh) 0.0417 (abayes) 
     *aster_h_t16 : 0.0473 (sh) 0.0403 (abayes) 
  *wtreeqmc_wh_n2 : 0.0431 (sh) 0.0382 (abayes) 
            *caml :             0.0552 (none) 
   wtreeqmc_wh_n1 : 0.0427 (sh) 0.0369 (abayes) 
   wtreeqmc_wh_n0 : 0.0427 (sh) 0.0369 (abayes) 
   wtreeqmc_ws_n2 : 0.0427 (sh) 0.0369 (abayes) 
   wtreeqmc_wn_n2 :             0.0469 (none) 

Model : 200 ntax, 2000000 hght, 1e-06 rate, 1000 genes
         asteroid :             0.9987 (none) 
        wastrid_s : 0.0304 (sh) 0.0260 (abayes) 
     *aster_h_t16 : 0.0303 (sh) 0.0241 (abayes) 
  *wtreeqmc_wh_n2 : 0.0262 (sh) 0.0227 (abayes) 
            *caml :             0.0278 (none) 
   wtreeqmc_wh_n1 : 0.0279 (sh) 0.0236 (abayes) 
   wtreeqmc_wh_n0 : 0.0279 (sh) 0.0236 (abayes) 
   wtreeqmc_ws_n2 : 0.0279 (sh) 0.0236 (abayes) 
   wtreeqmc_wn_n2 :             0.0310 (none) 

Model : 200 ntax, 500000 hght, 1e-07 rate, 50 genes
         asteroid :             0.9986 (none) 
        wastrid_s : 0.2154 (sh) 0.2063 (abayes) 
     *aster_h_t16 : 0.2019 (sh) 0.1898 (abayes) 
  *wtreeqmc_wh_n2 : 0.1956 (sh) 0.1858 (abayes) 
            *caml :             0.2816 (none) 
   wtreeqmc_wh_n1 : 0.1975 (sh) 0.1833 (abayes) 
   wtreeqmc_wh_n0 : 0.1975 (sh) 0.1833 (abayes) 
   wtreeqmc_ws_n2 : 0.1975 (sh) 0.1833 (abayes) 
   wtreeqmc_wn_n2 :             0.2042 (none) 

Model : 200 ntax, 500000 hght, 1e-07 rate, 200 genes
         asteroid :             0.9984 (none) 
        wastrid_s : 0.1133 (sh) 0.1056 (abayes) 
     *aster_h_t16 : 0.1016 (sh) 0.0913 (abayes) 
  *wtreeqmc_wh_n2 : 0.0981 (sh) 0.0919 (abayes) 
            *caml :             0.1610 (none) 
   wtreeqmc_wh_n1 : 0.0988 (sh) 0.0915 (abayes) 
   wtreeqmc_wh_n0 : 0.0988 (sh) 0.0915 (abayes) 
   wtreeqmc_ws_n2 : 0.0988 (sh) 0.0915 (abayes) 
   wtreeqmc_wn_n2 :             0.1076 (none) 

Model : 200 ntax, 500000 hght, 1e-07 rate, 1000 genes
         asteroid :             0.9985 (none) 
        wastrid_s : 0.0567 (sh) 0.0494 (abayes) 
     *aster_h_t16 : 0.0554 (sh) 0.0473 (abayes) 
  *wtreeqmc_wh_n2 : 0.0554 (sh) 0.0487 (abayes) 
            *caml :             0.0799 (none) 
   wtreeqmc_wh_n1 : 0.0562 (sh) 0.0492 (abayes) 
   wtreeqmc_wh_n0 : 0.0562 (sh) 0.0492 (abayes) 
   wtreeqmc_ws_n2 : 0.0562 (sh) 0.0492 (abayes) 
   wtreeqmc_wn_n2 :             0.0609 (none) 

Model : 200 ntax, 500000 hght, 1e-06 rate, 50 genes
         asteroid :             0.9991 (none) 
        wastrid_s : 0.2092 (sh) 0.2022 (abayes) 
     *aster_h_t16 : 0.1880 (sh) 0.1818 (abayes) 
  *wtreeqmc_wh_n2 : 0.1785 (sh) 0.1743 (abayes) 
            *caml :             0.2795 (none) 
   wtreeqmc_wh_n1 : 0.1789 (sh) 0.1706 (abayes) 
   wtreeqmc_wh_n0 : 0.1789 (sh) 0.1706 (abayes) 
   wtreeqmc_ws_n2 : 0.1789 (sh) 0.1706 (abayes) 
   wtreeqmc_wn_n2 :             0.1839 (none) 

Model : 200 ntax, 500000 hght, 1e-06 rate, 200 genes
         asteroid :             0.9992 (none) 
        wastrid_s : 0.1094 (sh) 0.1013 (abayes) 
     *aster_h_t16 : 0.0960 (sh) 0.0894 (abayes) 
  *wtreeqmc_wh_n2 : 0.0916 (sh) 0.0830 (abayes) 
            *caml :             0.1625 (none) 
   wtreeqmc_wh_n1 : 0.0888 (sh) 0.0811 (abayes) 
   wtreeqmc_wh_n0 : 0.0888 (sh) 0.0811 (abayes) 
   wtreeqmc_ws_n2 : 0.0888 (sh) 0.0811 (abayes) 
   wtreeqmc_wn_n2 :             0.0965 (none) 

Model : 200 ntax, 500000 hght, 1e-06 rate, 1000 genes
         asteroid :             0.9992 (none) 
        wastrid_s : 0.0512 (sh) 0.0471 (abayes) 
     *aster_h_t16 : 0.0471 (sh) 0.0400 (abayes) 
  *wtreeqmc_wh_n2 : 0.0436 (sh) 0.0377 (abayes) 
            *caml :             0.0795 (none) 
   wtreeqmc_wh_n1 : 0.0459 (sh) 0.0404 (abayes) 
   wtreeqmc_wh_n0 : 0.0459 (sh) 0.0404 (abayes) 
   wtreeqmc_ws_n2 : 0.0459 (sh) 0.0404 (abayes) 
   wtreeqmc_wn_n2 :             0.0493 (none) 
"""