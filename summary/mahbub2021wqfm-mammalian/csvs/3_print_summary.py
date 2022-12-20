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

for do in ["ils", "ngen", "slen"]: 
    if do == "ils":
        scals =["noscale", "scale2d", "scale2u"]
        ngens = ["200g"]
        nbpss = ["500b"]
    elif do == "ngen":
        scals = ["noscale"]
        ngens = ["25g", "50g", "100g", "400g", "800g"]
        nbpss = ["500b"]
    elif do == "slen":
        scals = ["noscale"]
        ngens = ["200g"]
        nbpss = ["250b", "1000b", "1500b", "true"]

    repls = [str("R%d" % x) for x in range(1, 21)]

    for scal in scals:
        for ngen in ngens:
            for nbps in nbpss:
                print("Model : %s, %s ngen, %s nbps" % (scal, ngen, nbps))
                xdf = df[(df["SCAL"] == scal) &
                         (df["NGEN"] == ngen) &
                         (df["NBPS"] == nbps)]

                for mthd, name in zip(mthds, names):
                    ydf = xdf[(xdf["MTHD"] == mthd)]
                    data = ydf.SERF.values

                    if data.size != 20:
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
Model : noscale, 200g ngen, 500b nbps
    treeqmc_n0_v1.0.0 : 0.0382 (err avg) 0.0294 (err med) 0.02 (mins)
    treeqmc_n1_v1.0.0 : 0.0382 (err avg) 0.0294 (err med) 0.02 (mins)
    treeqmc_n2_v1.0.0 : 0.0338 (err avg) 0.0294 (err med) 0.02 (mins)
      astral_3_v5.7.7 : 0.0382 (err avg) 0.0294 (err med) 0.05 (mins)
              fastral : 0.0382 (err avg) 0.0294 (err med) 0.09 (mins)
            wqfm_v1.3 : 0.0338 (err avg) 0.0294 (err med) 1.05 (mins)
            wqmc_v3.0 : 0.0382 (err avg) 0.0294 (err med) 0.85 (mins)

Model : scale2d, 200g ngen, 500b nbps
    treeqmc_n0_v1.0.0 : 0.0485 (err avg) 0.0441 (err med) 0.02 (mins)
    treeqmc_n1_v1.0.0 : 0.0485 (err avg) 0.0441 (err med) 0.02 (mins)
    treeqmc_n2_v1.0.0 : 0.0426 (err avg) 0.0294 (err med) 0.02 (mins)
      astral_3_v5.7.7 : 0.0485 (err avg) 0.0441 (err med) 0.10 (mins)
              fastral : 0.0485 (err avg) 0.0441 (err med) 0.07 (mins)
            wqfm_v1.3 : 0.0426 (err avg) 0.0294 (err med) 1.03 (mins)
            wqmc_v3.0 : 0.0559 (err avg) 0.0588 (err med) 0.87 (mins)

Model : scale2u, 200g ngen, 500b nbps
    treeqmc_n0_v1.0.0 : 0.0250 (err avg) 0.0294 (err med) 0.02 (mins)
    treeqmc_n1_v1.0.0 : 0.0250 (err avg) 0.0294 (err med) 0.02 (mins)
    treeqmc_n2_v1.0.0 : 0.0221 (err avg) 0.0294 (err med) 0.02 (mins)
      astral_3_v5.7.7 : 0.0250 (err avg) 0.0294 (err med) 0.03 (mins)
              fastral : 0.0250 (err avg) 0.0294 (err med) 0.07 (mins)
            wqfm_v1.3 : 0.0221 (err avg) 0.0294 (err med) 1.02 (mins)
            wqmc_v3.0 : 0.0250 (err avg) 0.0294 (err med) 0.87 (mins)

Model : noscale, 25g ngen, 500b nbps
    treeqmc_n0_v1.0.0 : 0.1132 (err avg) 0.1176 (err med) 0.00 (mins)
    treeqmc_n1_v1.0.0 : 0.1132 (err avg) 0.1176 (err med) 0.00 (mins)
    treeqmc_n2_v1.0.0 : 0.1118 (err avg) 0.1176 (err med) 0.00 (mins)
      astral_3_v5.7.7 : 0.1132 (err avg) 0.1176 (err med) 0.02 (mins)
              fastral : 0.1132 (err avg) 0.1176 (err med) 0.05 (mins)
            wqfm_v1.3 : 0.1103 (err avg) 0.1176 (err med) 0.28 (mins)
            wqmc_v3.0 : 0.1147 (err avg) 0.1176 (err med) 0.14 (mins)

Model : noscale, 50g ngen, 500b nbps
    treeqmc_n0_v1.0.0 : 0.0750 (err avg) 0.0882 (err med) 0.00 (mins)
    treeqmc_n1_v1.0.0 : 0.0809 (err avg) 0.0882 (err med) 0.00 (mins)
    treeqmc_n2_v1.0.0 : 0.0765 (err avg) 0.0588 (err med) 0.00 (mins)
      astral_3_v5.7.7 : 0.0750 (err avg) 0.0882 (err med) 0.02 (mins)
              fastral : 0.0750 (err avg) 0.0882 (err med) 0.05 (mins)
            wqfm_v1.3 : 0.0765 (err avg) 0.0588 (err med) 0.41 (mins)
            wqmc_v3.0 : 0.0779 (err avg) 0.0882 (err med) 0.24 (mins)

Model : noscale, 100g ngen, 500b nbps
    treeqmc_n0_v1.0.0 : 0.0706 (err avg) 0.0588 (err med) 0.01 (mins)
    treeqmc_n1_v1.0.0 : 0.0691 (err avg) 0.0588 (err med) 0.01 (mins)
    treeqmc_n2_v1.0.0 : 0.0647 (err avg) 0.0735 (err med) 0.01 (mins)
      astral_3_v5.7.7 : 0.0706 (err avg) 0.0588 (err med) 0.03 (mins)
              fastral : 0.0706 (err avg) 0.0588 (err med) 0.06 (mins)
            wqfm_v1.3 : 0.0647 (err avg) 0.0735 (err med) 0.63 (mins)
            wqmc_v3.0 : 0.0706 (err avg) 0.0588 (err med) 0.45 (mins)

Model : noscale, 400g ngen, 500b nbps
    treeqmc_n0_v1.0.0 : 0.0147 (err avg) 0.0000 (err med) 0.03 (mins)
    treeqmc_n1_v1.0.0 : 0.0147 (err avg) 0.0000 (err med) 0.04 (mins)
    treeqmc_n2_v1.0.0 : 0.0176 (err avg) 0.0147 (err med) 0.04 (mins)
      astral_3_v5.7.7 : 0.0147 (err avg) 0.0000 (err med) 0.08 (mins)
              fastral : 0.0147 (err avg) 0.0000 (err med) 0.10 (mins)
            wqfm_v1.3 : 0.0176 (err avg) 0.0147 (err med) 1.89 (mins)
            wqmc_v3.0 : 0.0147 (err avg) 0.0000 (err med) 1.70 (mins)

Model : noscale, 800g ngen, 500b nbps
    treeqmc_n0_v1.0.0 : 0.0074 (err avg) 0.0000 (err med) 0.07 (mins)
    treeqmc_n1_v1.0.0 : 0.0074 (err avg) 0.0000 (err med) 0.07 (mins)
    treeqmc_n2_v1.0.0 : 0.0059 (err avg) 0.0000 (err med) 0.07 (mins)
      astral_3_v5.7.7 : 0.0074 (err avg) 0.0000 (err med) 0.18 (mins)
              fastral : 0.0074 (err avg) 0.0000 (err med) 0.14 (mins)
            wqfm_v1.3 : 0.0059 (err avg) 0.0000 (err med) 3.53 (mins)
            wqmc_v3.0 : 0.0074 (err avg) 0.0000 (err med) 3.34 (mins)

Model : noscale, 200g ngen, 250b nbps
    treeqmc_n0_v1.0.0 : 0.0632 (err avg) 0.0588 (err med) 0.02 (mins)
    treeqmc_n1_v1.0.0 : 0.0647 (err avg) 0.0588 (err med) 0.02 (mins)
    treeqmc_n2_v1.0.0 : 0.0471 (err avg) 0.0294 (err med) 0.02 (mins)
      astral_3_v5.7.7 : 0.0632 (err avg) 0.0588 (err med) 0.07 (mins)
              fastral : 0.0632 (err avg) 0.0588 (err med) 0.07 (mins)
            wqfm_v1.3 : 0.0471 (err avg) 0.0294 (err med) 1.05 (mins)
            wqmc_v3.0 : 0.0691 (err avg) 0.0588 (err med) 0.86 (mins)

Model : noscale, 200g ngen, 1000b nbps
    treeqmc_n0_v1.0.0 : 0.0221 (err avg) 0.0294 (err med) 0.02 (mins)
    treeqmc_n1_v1.0.0 : 0.0206 (err avg) 0.0294 (err med) 0.02 (mins)
    treeqmc_n2_v1.0.0 : 0.0191 (err avg) 0.0294 (err med) 0.02 (mins)
      astral_3_v5.7.7 : 0.0221 (err avg) 0.0294 (err med) 0.04 (mins)
              fastral : 0.0221 (err avg) 0.0294 (err med) 0.07 (mins)
            wqfm_v1.3 : 0.0191 (err avg) 0.0294 (err med) 1.04 (mins)
            wqmc_v3.0 : 0.0221 (err avg) 0.0294 (err med) 0.86 (mins)

Model : noscale, 200g ngen, 1500b nbps
    treeqmc_n0_v1.0.0 : 0.0221 (err avg) 0.0294 (err med) 0.02 (mins)
    treeqmc_n1_v1.0.0 : 0.0162 (err avg) 0.0147 (err med) 0.02 (mins)
    treeqmc_n2_v1.0.0 : 0.0147 (err avg) 0.0000 (err med) 0.02 (mins)
      astral_3_v5.7.7 : 0.0176 (err avg) 0.0294 (err med) 0.04 (mins)
              fastral : 0.0176 (err avg) 0.0294 (err med) 0.07 (mins)
            wqfm_v1.3 : 0.0147 (err avg) 0.0000 (err med) 1.03 (mins)
            wqmc_v3.0 : 0.0221 (err avg) 0.0294 (err med) 0.86 (mins)

Model : noscale, 200g ngen, true nbps
    treeqmc_n0_v1.0.0 : 0.0132 (err avg) 0.0000 (err med) 0.02 (mins)
    treeqmc_n1_v1.0.0 : 0.0132 (err avg) 0.0000 (err med) 0.02 (mins)
    treeqmc_n2_v1.0.0 : 0.0103 (err avg) 0.0000 (err med) 0.02 (mins)
      astral_3_v5.7.7 : 0.0132 (err avg) 0.0000 (err med) 0.03 (mins)
              fastral : 0.0132 (err avg) 0.0000 (err med) 0.07 (mins)
            wqfm_v1.3 : 0.0103 (err avg) 0.0000 (err med) 1.04 (mins)
            wqmc_v3.0 : 0.0147 (err avg) 0.0000 (err med) 0.87 (mins)
"""