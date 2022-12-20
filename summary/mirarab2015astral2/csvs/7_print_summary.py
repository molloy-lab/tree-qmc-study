import pandas
import numpy
import sys

df = pandas.read_csv("data-all-error-and-timings-ex.csv", keep_default_na=False)

mthds = ["treeqmc_n0_v1.0.0", 
         "treeqmc_n1_v1.0.0", 
         "treeqmc_n2_v1.0.0",
         "astral_3_v5.7.7", 
         "fastral",
         "wqfm_v1.3",
         "wqmc_v3.0"]
names = ["treeqmc_n0_v1.0.0", 
         "treeqmc_n1_v1.0.0", 
         "treeqmc_n2_v1.0.0",
         "  astral_3_v5.7.7", 
         "          fastral",
         "        wqfm_v1.3",
         "        wqmc_v3.0"]


for do in ["ntax", "ils"]: 
    if do == "ntax":
        ntaxs = [10, 50, 100, 500, 1000]
        strhts = [2000000]
        srates = [0.000001]
        
    elif do == "ils":
        ntaxs = [200]
        strhts = [10000000, 2000000, 500000]
        srates = [0.0000001, 0.000001]
    gtre = "estimatedgenetre"  # "truegenetrees"
    ngen = 1000

    for ntax in ntaxs:
        for strht in strhts:
            for srate in srates:
                print("Model : %s ntax, %s strht, %s srate" % (ntax, strht, srate))
                xdf = df[(df["NTAX"] == ntax) &
                         (df["STRHT"] == strht) &
                         (df["SRATE"] == srate) &
                         (df["GTRE"] == gtre) &
                         (df["NGEN"] == ngen)]

                for mthd, name in zip(mthds, names):
                    ydf = xdf[(xdf["MTHD"] == mthd)]

                    data = ydf.SERF.values
                    rfavg = numpy.mean(data)
                    rfavg = numpy.round(numpy.round(rfavg, 5), 4)

                    rfmed = numpy.median(data)
                    rfmed = numpy.round(numpy.round(rfmed, 5), 4)

                    data = ydf.SECS.values / 60
                    mins = numpy.mean(data)
                    mins = numpy.round(numpy.round(mins, 3), 2)

                    print("    %s : %1.4f (err avg) %1.4f (err med) %1.2f (mins)" % (name, rfavg, rfmed, mins))

                print()

"""
Model : 10 ntax, 2000000 strht, 1e-06 srate
    treeqmc_n0_v1.0.0 : 0.0153 (err avg) 0.0000 (err med) 0.01 (mins)
    treeqmc_n1_v1.0.0 : 0.0153 (err avg) 0.0000 (err med) 0.01 (mins)
    treeqmc_n2_v1.0.0 : 0.0153 (err avg) 0.0000 (err med) 0.01 (mins)
      astral_3_v5.7.7 : 0.0153 (err avg) 0.0000 (err med) 0.03 (mins)
              fastral : 0.0153 (err avg) 0.0000 (err med) 0.12 (mins)
            wqfm_v1.3 : 0.0153 (err avg) 0.0000 (err med) 0.23 (mins)
            wqmc_v3.0 : 0.0153 (err avg) 0.0000 (err med) 0.21 (mins)

Model : 50 ntax, 2000000 strht, 1e-06 srate
    treeqmc_n0_v1.0.0 : 0.0278 (err avg) 0.0208 (err med) 0.17 (mins)
    treeqmc_n1_v1.0.0 : 0.0278 (err avg) 0.0208 (err med) 0.17 (mins)
    treeqmc_n2_v1.0.0 : 0.0260 (err avg) 0.0208 (err med) 0.16 (mins)
      astral_3_v5.7.7 : 0.0269 (err avg) 0.0208 (err med) 0.60 (mins)
              fastral : 0.0269 (err avg) 0.0208 (err med) 0.29 (mins)
            wqfm_v1.3 : 0.0265 (err avg) 0.0208 (err med) 16.79 (mins)
            wqmc_v3.0 : 0.0278 (err avg) 0.0208 (err med) 15.83 (mins)

Model : 100 ntax, 2000000 strht, 1e-06 srate
    treeqmc_n0_v1.0.0 : 0.0334 (err avg) 0.0204 (err med) 0.60 (mins)
    treeqmc_n1_v1.0.0 : 0.0319 (err avg) 0.0204 (err med) 0.63 (mins)
    treeqmc_n2_v1.0.0 : 0.0306 (err avg) 0.0204 (err med) 0.62 (mins)
      astral_3_v5.7.7 : 0.0298 (err avg) 0.0204 (err med) 3.05 (mins)
              fastral : 0.0310 (err avg) 0.0204 (err med) 0.67 (mins)
            wqfm_v1.3 : 0.0296 (err avg) 0.0204 (err med) 290.16 (mins)
            wqmc_v3.0 : 0.0370 (err avg) 0.0204 (err med) 256.75 (mins)

Model : 500 ntax, 2000000 strht, 1e-06 srate
    treeqmc_n0_v1.0.0 : 0.0392 (err avg) 0.0311 (err med) 14.83 (mins)
    treeqmc_n1_v1.0.0 : 0.0316 (err avg) 0.0251 (err med) 15.82 (mins)
    treeqmc_n2_v1.0.0 : 0.0280 (err avg) 0.0231 (err med) 15.48 (mins)
      astral_3_v5.7.7 : 0.0334 (err avg) 0.0251 (err med) 72.08 (mins)
              fastral : 0.0328 (err avg) 0.0251 (err med) 5.59 (mins)
            wqfm_v1.3 : nan (err avg) nan (err med) nan (mins)
            wqmc_v3.0 : nan (err avg) nan (err med) nan (mins)

Model : 1000 ntax, 2000000 strht, 1e-06 srate
    treeqmc_n0_v1.0.0 : 0.0518 (err avg) 0.0291 (err med) 59.15 (mins)
    treeqmc_n1_v1.0.0 : 0.0407 (err avg) 0.0240 (err med) 62.03 (mins)
    treeqmc_n2_v1.0.0 : 0.0360 (err avg) 0.0200 (err med) 63.51 (mins)
      astral_3_v5.7.7 : 0.0407 (err avg) 0.0240 (err med) 317.77 (mins)
              fastral : 0.0400 (err avg) 0.0240 (err med) 32.37 (mins)
            wqfm_v1.3 : nan (err avg) nan (err med) nan (mins)
            wqmc_v3.0 : nan (err avg) nan (err med) nan (mins)

Model : 200 ntax, 10000000 strht, 1e-07 srate
    treeqmc_n0_v1.0.0 : 0.0562 (err avg) 0.0328 (err med) 2.22 (mins)
    treeqmc_n1_v1.0.0 : 0.0486 (err avg) 0.0303 (err med) 2.34 (mins)
    treeqmc_n2_v1.0.0 : 0.0450 (err avg) 0.0253 (err med) 2.34 (mins)
      astral_3_v5.7.7 : 0.0496 (err avg) 0.0354 (err med) 6.56 (mins)
              fastral : 0.0506 (err avg) 0.0354 (err med) 1.17 (mins)
            wqfm_v1.3 : nan (err avg) nan (err med) nan (mins)
            wqmc_v3.0 : nan (err avg) nan (err med) nan (mins)

Model : 200 ntax, 10000000 strht, 1e-06 srate
    treeqmc_n0_v1.0.0 : 0.0206 (err avg) 0.0101 (err med) 2.24 (mins)
    treeqmc_n1_v1.0.0 : 0.0203 (err avg) 0.0101 (err med) 2.31 (mins)
    treeqmc_n2_v1.0.0 : 0.0183 (err avg) 0.0101 (err med) 2.31 (mins)
      astral_3_v5.7.7 : 0.0204 (err avg) 0.0101 (err med) 4.64 (mins)
              fastral : 0.0205 (err avg) 0.0101 (err med) 1.14 (mins)
            wqfm_v1.3 : nan (err avg) nan (err med) nan (mins)
            wqmc_v3.0 : nan (err avg) nan (err med) nan (mins)

Model : 200 ntax, 2000000 strht, 1e-07 srate
    treeqmc_n0_v1.0.0 : 0.0501 (err avg) 0.0303 (err med) 2.37 (mins)
    treeqmc_n1_v1.0.0 : 0.0439 (err avg) 0.0303 (err med) 2.54 (mins)
    treeqmc_n2_v1.0.0 : 0.0396 (err avg) 0.0253 (err med) 2.53 (mins)
      astral_3_v5.7.7 : 0.0446 (err avg) 0.0278 (err med) 17.11 (mins)
              fastral : 0.0431 (err avg) 0.0278 (err med) 1.24 (mins)
            wqfm_v1.3 : nan (err avg) nan (err med) nan (mins)
            wqmc_v3.0 : nan (err avg) nan (err med) nan (mins)

Model : 200 ntax, 2000000 strht, 1e-06 srate
    treeqmc_n0_v1.0.0 : 0.0364 (err avg) 0.0253 (err med) 2.35 (mins)
    treeqmc_n1_v1.0.0 : 0.0334 (err avg) 0.0253 (err med) 2.48 (mins)
    treeqmc_n2_v1.0.0 : 0.0315 (err avg) 0.0202 (err med) 2.47 (mins)
      astral_3_v5.7.7 : 0.0339 (err avg) 0.0253 (err med) 12.05 (mins)
              fastral : 0.0337 (err avg) 0.0253 (err med) 1.24 (mins)
            wqfm_v1.3 : nan (err avg) nan (err med) nan (mins)
            wqmc_v3.0 : nan (err avg) nan (err med) nan (mins)

Model : 200 ntax, 500000 strht, 1e-07 srate
    treeqmc_n0_v1.0.0 : 0.0923 (err avg) 0.0758 (err med) 2.66 (mins)
    treeqmc_n1_v1.0.0 : 0.0668 (err avg) 0.0505 (err med) 2.95 (mins)
    treeqmc_n2_v1.0.0 : 0.0628 (err avg) 0.0455 (err med) 2.96 (mins)
      astral_3_v5.7.7 : 0.0654 (err avg) 0.0505 (err med) 77.36 (mins)
              fastral : 0.0616 (err avg) 0.0505 (err med) 1.58 (mins)
            wqfm_v1.3 : nan (err avg) nan (err med) nan (mins)
            wqmc_v3.0 : nan (err avg) nan (err med) nan (mins)

Model : 200 ntax, 500000 strht, 1e-06 srate
    treeqmc_n0_v1.0.0 : 0.0642 (err avg) 0.0556 (err med) 2.59 (mins)
    treeqmc_n1_v1.0.0 : 0.0522 (err avg) 0.0455 (err med) 2.84 (mins)
    treeqmc_n2_v1.0.0 : 0.0498 (err avg) 0.0404 (err med) 2.78 (mins)
      astral_3_v5.7.7 : 0.0536 (err avg) 0.0505 (err med) 72.71 (mins)
              fastral : 0.0521 (err avg) 0.0455 (err med) 1.50 (mins)
            wqfm_v1.3 : nan (err avg) nan (err med) nan (mins)
            wqmc_v3.0 : nan (err avg) nan (err med) nan (mins)
"""