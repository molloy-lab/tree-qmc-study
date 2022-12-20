import pandas
import numpy
import sys

df = pandas.read_csv("data-all-error-and-timings.csv", keep_default_na=False)

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

for do in ["ils", "ngen"]:
    if do == "ils":
        scals =["0.5X", "1X", "2X"]
        ngens = [1000]
        nbpss = ["true", "500"]
    elif do == "ngen":
        scals = ["1X"]
        ngens = [50, 100, 200, 500]
        nbpss = ["true", "500"]

    for scal in scals:
        for ngen in ngens:
            for nbps in nbpss:
                nrepl = 20

                print("Model : %s, %s ngen, %s nbps" % (scal, ngen, nbps))
                xdf = df[(df["SCAL"] == scal) &
                         (df["NGEN"] == ngen) &
                         (df["NBPS"] == nbps)]

                for mthd, name in zip(mthds, names):
                    ydf = xdf[(xdf["MTHD"] == mthd)]
                    data = ydf.SERF.values

                    if data.size != nrepl:
                        sys.exit("ERROR - wrong number of replicates!")

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
Model : 0.5X, 1000 ngen, true nbps
    treeqmc_n0_v1.0.0 : 0.0478 (err avg) 0.0444 (err med) 0.18 (mins)
    treeqmc_n1_v1.0.0 : 0.0411 (err avg) 0.0444 (err med) 0.20 (mins)
    treeqmc_n2_v1.0.0 : 0.0444 (err avg) 0.0444 (err med) 0.19 (mins)
      astral_3_v5.7.7 : 0.0478 (err avg) 0.0444 (err med) 3.16 (mins)
              fastral : 0.0400 (err avg) 0.0444 (err med) 0.30 (mins)
            wqfm_v1.3 : 0.0433 (err avg) 0.0444 (err med) 12.39 (mins)
            wqmc_v3.0 : 0.0944 (err avg) 0.0889 (err med) 11.66 (mins)

Model : 0.5X, 1000 ngen, 500 nbps
    treeqmc_n0_v1.0.0 : 0.1200 (err avg) 0.1111 (err med) 0.19 (mins)
    treeqmc_n1_v1.0.0 : 0.0978 (err avg) 0.0889 (err med) 0.20 (mins)
    treeqmc_n2_v1.0.0 : 0.0844 (err avg) 0.0778 (err med) 0.20 (mins)
      astral_3_v5.7.7 : 0.1311 (err avg) 0.1333 (err med) 4.32 (mins)
              fastral : 0.1011 (err avg) 0.1000 (err med) 0.35 (mins)
            wqfm_v1.3 : 0.1011 (err avg) 0.0889 (err med) 12.20 (mins)
            wqmc_v3.0 : 0.2222 (err avg) 0.2222 (err med) 11.48 (mins)

Model : 1X, 1000 ngen, true nbps
    treeqmc_n0_v1.0.0 : 0.0367 (err avg) 0.0222 (err med) 0.18 (mins)
    treeqmc_n1_v1.0.0 : 0.0300 (err avg) 0.0222 (err med) 0.19 (mins)
    treeqmc_n2_v1.0.0 : 0.0344 (err avg) 0.0444 (err med) 0.19 (mins)
      astral_3_v5.7.7 : 0.0333 (err avg) 0.0222 (err med) 2.11 (mins)
              fastral : 0.0278 (err avg) 0.0222 (err med) 0.32 (mins)
            wqfm_v1.3 : 0.0289 (err avg) 0.0222 (err med) 12.81 (mins)
            wqmc_v3.0 : 0.0600 (err avg) 0.0667 (err med) 12.15 (mins)

Model : 1X, 1000 ngen, 500 nbps
    treeqmc_n0_v1.0.0 : 0.0956 (err avg) 0.1000 (err med) 0.17 (mins)
    treeqmc_n1_v1.0.0 : 0.0567 (err avg) 0.0556 (err med) 0.18 (mins)
    treeqmc_n2_v1.0.0 : 0.0533 (err avg) 0.0556 (err med) 0.18 (mins)
      astral_3_v5.7.7 : 0.0778 (err avg) 0.0667 (err med) 4.06 (mins)
              fastral : 0.0689 (err avg) 0.0778 (err med) 0.29 (mins)
            wqfm_v1.3 : 0.0544 (err avg) 0.0444 (err med) 12.36 (mins)
            wqmc_v3.0 : 0.1333 (err avg) 0.1333 (err med) 11.60 (mins)

Model : 2X, 1000 ngen, true nbps
    treeqmc_n0_v1.0.0 : 0.0189 (err avg) 0.0222 (err med) 0.17 (mins)
    treeqmc_n1_v1.0.0 : 0.0222 (err avg) 0.0222 (err med) 0.18 (mins)
    treeqmc_n2_v1.0.0 : 0.0200 (err avg) 0.0222 (err med) 0.18 (mins)
      astral_3_v5.7.7 : 0.0211 (err avg) 0.0222 (err med) 1.27 (mins)
              fastral : 0.0211 (err avg) 0.0222 (err med) 0.28 (mins)
            wqfm_v1.3 : 0.0200 (err avg) 0.0222 (err med) 12.12 (mins)
            wqmc_v3.0 : 0.0333 (err avg) 0.0222 (err med) 11.52 (mins)

Model : 2X, 1000 ngen, 500 nbps
    treeqmc_n0_v1.0.0 : 0.0622 (err avg) 0.0667 (err med) 0.16 (mins)
    treeqmc_n1_v1.0.0 : 0.0422 (err avg) 0.0444 (err med) 0.17 (mins)
    treeqmc_n2_v1.0.0 : 0.0400 (err avg) 0.0444 (err med) 0.17 (mins)
      astral_3_v5.7.7 : 0.0433 (err avg) 0.0444 (err med) 3.51 (mins)
              fastral : 0.0444 (err avg) 0.0444 (err med) 0.19 (mins)
            wqfm_v1.3 : 0.0400 (err avg) 0.0444 (err med) 12.50 (mins)
            wqmc_v3.0 : 0.0922 (err avg) 0.0889 (err med) 11.71 (mins)

Model : 1X, 50 ngen, true nbps
    treeqmc_n0_v1.0.0 : 0.1900 (err avg) 0.2000 (err med) 0.01 (mins)
    treeqmc_n1_v1.0.0 : 0.1544 (err avg) 0.1556 (err med) 0.01 (mins)
    treeqmc_n2_v1.0.0 : 0.1544 (err avg) 0.1556 (err med) 0.01 (mins)
      astral_3_v5.7.7 : 0.1544 (err avg) 0.1556 (err med) 0.05 (mins)
              fastral : 0.1700 (err avg) 0.1778 (err med) 0.06 (mins)
            wqfm_v1.3 : 0.1400 (err avg) 0.1333 (err med) 1.26 (mins)
            wqmc_v3.0 : 0.1900 (err avg) 0.2000 (err med) 0.64 (mins)

Model : 1X, 50 ngen, 500 nbps
    treeqmc_n0_v1.0.0 : 0.3200 (err avg) 0.3111 (err med) 0.01 (mins)
    treeqmc_n1_v1.0.0 : 0.2844 (err avg) 0.2889 (err med) 0.01 (mins)
    treeqmc_n2_v1.0.0 : 0.2756 (err avg) 0.2667 (err med) 0.01 (mins)
      astral_3_v5.7.7 : 0.2833 (err avg) 0.2889 (err med) 0.07 (mins)
              fastral : 0.2933 (err avg) 0.2889 (err med) 0.07 (mins)
            wqfm_v1.3 : 0.2678 (err avg) 0.2667 (err med) 1.37 (mins)
            wqmc_v3.0 : 0.3144 (err avg) 0.3111 (err med) 0.64 (mins)

Model : 1X, 100 ngen, true nbps
    treeqmc_n0_v1.0.0 : 0.1356 (err avg) 0.1333 (err med) 0.02 (mins)
    treeqmc_n1_v1.0.0 : 0.1122 (err avg) 0.1111 (err med) 0.02 (mins)
    treeqmc_n2_v1.0.0 : 0.1078 (err avg) 0.1111 (err med) 0.02 (mins)
      astral_3_v5.7.7 : 0.1211 (err avg) 0.1111 (err med) 0.07 (mins)
              fastral : 0.1233 (err avg) 0.1111 (err med) 0.08 (mins)
            wqfm_v1.3 : 0.1033 (err avg) 0.1111 (err med) 1.82 (mins)
            wqmc_v3.0 : 0.1467 (err avg) 0.1333 (err med) 1.19 (mins)

Model : 1X, 100 ngen, 500 nbps
    treeqmc_n0_v1.0.0 : 0.2711 (err avg) 0.2667 (err med) 0.02 (mins)
    treeqmc_n1_v1.0.0 : 0.2311 (err avg) 0.2444 (err med) 0.02 (mins)
    treeqmc_n2_v1.0.0 : 0.2278 (err avg) 0.2333 (err med) 0.02 (mins)
      astral_3_v5.7.7 : 0.2522 (err avg) 0.2444 (err med) 0.12 (mins)
              fastral : 0.2567 (err avg) 0.2556 (err med) 0.08 (mins)
            wqfm_v1.3 : 0.2289 (err avg) 0.2444 (err med) 1.97 (mins)
            wqmc_v3.0 : 0.2833 (err avg) 0.2778 (err med) 1.25 (mins)

Model : 1X, 200 ngen, true nbps
    treeqmc_n0_v1.0.0 : 0.0833 (err avg) 0.0889 (err med) 0.03 (mins)
    treeqmc_n1_v1.0.0 : 0.0789 (err avg) 0.0667 (err med) 0.03 (mins)
    treeqmc_n2_v1.0.0 : 0.0700 (err avg) 0.0667 (err med) 0.03 (mins)
      astral_3_v5.7.7 : 0.0778 (err avg) 0.0889 (err med) 0.13 (mins)
              fastral : 0.0744 (err avg) 0.0778 (err med) 0.10 (mins)
            wqfm_v1.3 : 0.0678 (err avg) 0.0667 (err med) 2.96 (mins)
            wqmc_v3.0 : 0.1044 (err avg) 0.1111 (err med) 2.33 (mins)

Model : 1X, 200 ngen, 500 nbps
    treeqmc_n0_v1.0.0 : 0.2178 (err avg) 0.2222 (err med) 0.03 (mins)
    treeqmc_n1_v1.0.0 : 0.1722 (err avg) 0.1556 (err med) 0.03 (mins)
    treeqmc_n2_v1.0.0 : 0.1722 (err avg) 0.1778 (err med) 0.03 (mins)
      astral_3_v5.7.7 : 0.1944 (err avg) 0.1778 (err med) 0.24 (mins)
              fastral : 0.1789 (err avg) 0.1667 (err med) 0.11 (mins)
            wqfm_v1.3 : 0.1778 (err avg) 0.1667 (err med) 3.11 (mins)
            wqmc_v3.0 : 0.2322 (err avg) 0.2333 (err med) 2.39 (mins)

Model : 1X, 500 ngen, true nbps
    treeqmc_n0_v1.0.0 : 0.0611 (err avg) 0.0667 (err med) 0.08 (mins)
    treeqmc_n1_v1.0.0 : 0.0444 (err avg) 0.0444 (err med) 0.09 (mins)
    treeqmc_n2_v1.0.0 : 0.0389 (err avg) 0.0444 (err med) 0.09 (mins)
      astral_3_v5.7.7 : 0.0489 (err avg) 0.0444 (err med) 0.52 (mins)
              fastral : 0.0489 (err avg) 0.0556 (err med) 0.18 (mins)
            wqfm_v1.3 : 0.0456 (err avg) 0.0444 (err med) 6.44 (mins)
            wqmc_v3.0 : 0.0767 (err avg) 0.0667 (err med) 5.80 (mins)

Model : 1X, 500 ngen, 500 nbps
    treeqmc_n0_v1.0.0 : 0.1389 (err avg) 0.1333 (err med) 0.08 (mins)
    treeqmc_n1_v1.0.0 : 0.0967 (err avg) 0.0889 (err med) 0.09 (mins)
    treeqmc_n2_v1.0.0 : 0.0900 (err avg) 0.0889 (err med) 0.09 (mins)
      astral_3_v5.7.7 : 0.1233 (err avg) 0.1333 (err med) 1.02 (mins)
              fastral : 0.1133 (err avg) 0.1111 (err med) 0.18 (mins)
            wqfm_v1.3 : 0.1100 (err avg) 0.1111 (err med) 6.58 (mins)
            wqmc_v3.0 : 0.1800 (err avg) 0.1778 (err med) 5.84 (mins)
"""