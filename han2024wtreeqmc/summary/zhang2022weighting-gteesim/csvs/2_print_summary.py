import pandas
import numpy
import sys

sys.exit("DONE RUNNING")

df = pandas.read_csv("data-all-error.csv", na_values='NA', keep_default_na=False)

mthds = ["wastrid_s",
         "aster_h_t16",
         "wtreeqmc_wh_n2",
         "wtreeqmc_wh_n1",
         "wtreeqmc_wh_n0",
         "wtreeqmc_ws_n2",
         "wtreeqmc_wn_n2"]

names = ["       wastrid_s",
         "     aster_h_t16",
         "  wtreeqmc_wh_n2",
         "  wtreeqmc_wh_n1",
         "  wtreeqmc_wh_n0",
         "  wtreeqmc_ws_n2",
         "  wtreeqmc_wn_n2"]

nbpss = [200, 400, 800, 1600]
ngens = [50, 200, 500, 1000]
supps = ["bs", "abayes"]

for nbps in nbpss:
    for ngen in ngens:
        for supp in supps:
            print("%d bp, %d genes, %s support" % (nbps, ngen, supp))

            xdf = df[(df["NBPS"] == nbps) &
                     (df["NGEN"] == ngen) &
                     (df["SUPP"] == supp)]

            for name, mthd in zip(names, mthds):
                ydf = xdf[(xdf["MTHD"] == mthd)]
                data = ydf.SERF.values
                if data.size != 50:
                    sys.stdout.write("ERROR - wrong number of replicates!\n")

                nrepl = numpy.sum(~numpy.isnan(data))

                fnmed = numpy.mean(ydf.SEFN.values)
                fnmed = numpy.round(numpy.round(fnmed, 3), 2)

                rfmed = numpy.mean(ydf.SERF.values)
                rfmed = numpy.round(numpy.round(rfmed, 3), 2)

                print("    %s : %1.2f (fn avg) %1.2f (rf avg) - %d repls" % (name, fnmed, rfmed, nrepl))

            print()

"""
200 bp, 50 genes, bs support
           wastrid_s : 15.98 (fn avg) 0.16 (rf avg) - 50 repls
         aster_h_t16 : 15.22 (fn avg) 0.16 (rf avg) - 50 repls
      wtreeqmc_wh_n2 : 15.08 (fn avg) 0.15 (rf avg) - 50 repls
      wtreeqmc_wh_n1 : 15.30 (fn avg) 0.16 (rf avg) - 50 repls
      wtreeqmc_wh_n0 : 15.30 (fn avg) 0.16 (rf avg) - 50 repls
      wtreeqmc_ws_n2 : 15.30 (fn avg) 0.16 (rf avg) - 50 repls
      wtreeqmc_wn_n2 : 16.86 (fn avg) 0.17 (rf avg) - 50 repls

200 bp, 50 genes, abayes support
           wastrid_s : 14.74 (fn avg) 0.15 (rf avg) - 50 repls
         aster_h_t16 : 14.28 (fn avg) 0.15 (rf avg) - 50 repls
      wtreeqmc_wh_n2 : 13.86 (fn avg) 0.14 (rf avg) - 50 repls
      wtreeqmc_wh_n1 : 14.26 (fn avg) 0.15 (rf avg) - 50 repls
      wtreeqmc_wh_n0 : 14.26 (fn avg) 0.15 (rf avg) - 50 repls
      wtreeqmc_ws_n2 : 14.26 (fn avg) 0.15 (rf avg) - 50 repls
      wtreeqmc_wn_n2 : 17.04 (fn avg) 0.17 (rf avg) - 50 repls

200 bp, 200 genes, bs support
           wastrid_s : 9.52 (fn avg) 0.10 (rf avg) - 50 repls
         aster_h_t16 : 9.44 (fn avg) 0.10 (rf avg) - 50 repls
      wtreeqmc_wh_n2 : 9.30 (fn avg) 0.10 (rf avg) - 50 repls
      wtreeqmc_wh_n1 : 9.62 (fn avg) 0.10 (rf avg) - 50 repls
      wtreeqmc_wh_n0 : 9.62 (fn avg) 0.10 (rf avg) - 50 repls
      wtreeqmc_ws_n2 : 9.62 (fn avg) 0.10 (rf avg) - 50 repls
      wtreeqmc_wn_n2 : 10.72 (fn avg) 0.11 (rf avg) - 50 repls

200 bp, 200 genes, abayes support
           wastrid_s : 9.28 (fn avg) 0.10 (rf avg) - 50 repls
         aster_h_t16 : 9.04 (fn avg) 0.09 (rf avg) - 50 repls
      wtreeqmc_wh_n2 : 8.14 (fn avg) 0.08 (rf avg) - 50 repls
      wtreeqmc_wh_n1 : 9.24 (fn avg) 0.09 (rf avg) - 50 repls
      wtreeqmc_wh_n0 : 9.24 (fn avg) 0.09 (rf avg) - 50 repls
      wtreeqmc_ws_n2 : 9.24 (fn avg) 0.09 (rf avg) - 50 repls
      wtreeqmc_wn_n2 : 10.92 (fn avg) 0.11 (rf avg) - 50 repls

200 bp, 500 genes, bs support
           wastrid_s : 7.58 (fn avg) 0.08 (rf avg) - 50 repls
         aster_h_t16 : 7.38 (fn avg) 0.08 (rf avg) - 50 repls
      wtreeqmc_wh_n2 : 7.40 (fn avg) 0.08 (rf avg) - 50 repls
      wtreeqmc_wh_n1 : 7.78 (fn avg) 0.08 (rf avg) - 50 repls
      wtreeqmc_wh_n0 : 7.78 (fn avg) 0.08 (rf avg) - 50 repls
      wtreeqmc_ws_n2 : 7.78 (fn avg) 0.08 (rf avg) - 50 repls
      wtreeqmc_wn_n2 : 9.10 (fn avg) 0.09 (rf avg) - 50 repls

200 bp, 500 genes, abayes support
           wastrid_s : 7.38 (fn avg) 0.08 (rf avg) - 50 repls
         aster_h_t16 : 7.02 (fn avg) 0.07 (rf avg) - 50 repls
      wtreeqmc_wh_n2 : 6.78 (fn avg) 0.07 (rf avg) - 50 repls
      wtreeqmc_wh_n1 : 7.68 (fn avg) 0.08 (rf avg) - 50 repls
      wtreeqmc_wh_n0 : 7.68 (fn avg) 0.08 (rf avg) - 50 repls
      wtreeqmc_ws_n2 : 7.68 (fn avg) 0.08 (rf avg) - 50 repls
      wtreeqmc_wn_n2 : 9.06 (fn avg) 0.09 (rf avg) - 50 repls

200 bp, 1000 genes, bs support
           wastrid_s : 6.56 (fn avg) 0.07 (rf avg) - 50 repls
         aster_h_t16 : 6.40 (fn avg) 0.06 (rf avg) - 50 repls
      wtreeqmc_wh_n2 : 6.60 (fn avg) 0.07 (rf avg) - 50 repls
      wtreeqmc_wh_n1 : 6.90 (fn avg) 0.07 (rf avg) - 50 repls
      wtreeqmc_wh_n0 : 6.90 (fn avg) 0.07 (rf avg) - 50 repls
      wtreeqmc_ws_n2 : 6.90 (fn avg) 0.07 (rf avg) - 50 repls
      wtreeqmc_wn_n2 : 7.80 (fn avg) 0.08 (rf avg) - 50 repls

200 bp, 1000 genes, abayes support
           wastrid_s : 6.10 (fn avg) 0.06 (rf avg) - 50 repls
         aster_h_t16 : 6.18 (fn avg) 0.06 (rf avg) - 50 repls
      wtreeqmc_wh_n2 : 5.24 (fn avg) 0.05 (rf avg) - 50 repls
      wtreeqmc_wh_n1 : 5.80 (fn avg) 0.06 (rf avg) - 50 repls
      wtreeqmc_wh_n0 : 5.80 (fn avg) 0.06 (rf avg) - 50 repls
      wtreeqmc_ws_n2 : 5.80 (fn avg) 0.06 (rf avg) - 50 repls
      wtreeqmc_wn_n2 : 7.84 (fn avg) 0.08 (rf avg) - 50 repls

400 bp, 50 genes, bs support
           wastrid_s : 12.86 (fn avg) 0.13 (rf avg) - 50 repls
         aster_h_t16 : 12.30 (fn avg) 0.13 (rf avg) - 50 repls
      wtreeqmc_wh_n2 : 11.98 (fn avg) 0.12 (rf avg) - 50 repls
      wtreeqmc_wh_n1 : 12.10 (fn avg) 0.12 (rf avg) - 50 repls
      wtreeqmc_wh_n0 : 12.10 (fn avg) 0.12 (rf avg) - 50 repls
      wtreeqmc_ws_n2 : 12.10 (fn avg) 0.12 (rf avg) - 50 repls
      wtreeqmc_wn_n2 : 12.80 (fn avg) 0.13 (rf avg) - 50 repls

400 bp, 50 genes, abayes support
           wastrid_s : 12.14 (fn avg) 0.12 (rf avg) - 50 repls
         aster_h_t16 : 11.94 (fn avg) 0.12 (rf avg) - 50 repls
      wtreeqmc_wh_n2 : 10.90 (fn avg) 0.11 (rf avg) - 50 repls
      wtreeqmc_wh_n1 : 11.30 (fn avg) 0.12 (rf avg) - 50 repls
      wtreeqmc_wh_n0 : 11.30 (fn avg) 0.12 (rf avg) - 50 repls
      wtreeqmc_ws_n2 : 11.30 (fn avg) 0.12 (rf avg) - 50 repls
      wtreeqmc_wn_n2 : 12.86 (fn avg) 0.13 (rf avg) - 50 repls

400 bp, 200 genes, bs support
           wastrid_s : 7.80 (fn avg) 0.08 (rf avg) - 50 repls
         aster_h_t16 : 7.52 (fn avg) 0.08 (rf avg) - 50 repls
      wtreeqmc_wh_n2 : 7.14 (fn avg) 0.07 (rf avg) - 50 repls
      wtreeqmc_wh_n1 : 7.32 (fn avg) 0.08 (rf avg) - 50 repls
      wtreeqmc_wh_n0 : 7.32 (fn avg) 0.08 (rf avg) - 50 repls
      wtreeqmc_ws_n2 : 7.32 (fn avg) 0.08 (rf avg) - 50 repls
      wtreeqmc_wn_n2 : 8.36 (fn avg) 0.08 (rf avg) - 50 repls

400 bp, 200 genes, abayes support
           wastrid_s : 7.28 (fn avg) 0.07 (rf avg) - 50 repls
         aster_h_t16 : 7.08 (fn avg) 0.07 (rf avg) - 50 repls
      wtreeqmc_wh_n2 : 6.58 (fn avg) 0.07 (rf avg) - 50 repls
      wtreeqmc_wh_n1 : 7.04 (fn avg) 0.07 (rf avg) - 50 repls
      wtreeqmc_wh_n0 : 7.04 (fn avg) 0.07 (rf avg) - 50 repls
      wtreeqmc_ws_n2 : 7.04 (fn avg) 0.07 (rf avg) - 50 repls
      wtreeqmc_wn_n2 : 8.42 (fn avg) 0.09 (rf avg) - 50 repls

400 bp, 500 genes, bs support
           wastrid_s : 6.04 (fn avg) 0.06 (rf avg) - 50 repls
         aster_h_t16 : 5.84 (fn avg) 0.06 (rf avg) - 50 repls
      wtreeqmc_wh_n2 : 5.44 (fn avg) 0.06 (rf avg) - 50 repls
      wtreeqmc_wh_n1 : 5.88 (fn avg) 0.06 (rf avg) - 50 repls
      wtreeqmc_wh_n0 : 5.88 (fn avg) 0.06 (rf avg) - 50 repls
      wtreeqmc_ws_n2 : 5.88 (fn avg) 0.06 (rf avg) - 50 repls
      wtreeqmc_wn_n2 : 7.04 (fn avg) 0.07 (rf avg) - 50 repls

400 bp, 500 genes, abayes support
           wastrid_s : 5.96 (fn avg) 0.06 (rf avg) - 50 repls
         aster_h_t16 : 5.64 (fn avg) 0.06 (rf avg) - 50 repls
      wtreeqmc_wh_n2 : 5.74 (fn avg) 0.06 (rf avg) - 50 repls
      wtreeqmc_wh_n1 : 6.04 (fn avg) 0.06 (rf avg) - 50 repls
      wtreeqmc_wh_n0 : 6.04 (fn avg) 0.06 (rf avg) - 50 repls
      wtreeqmc_ws_n2 : 6.04 (fn avg) 0.06 (rf avg) - 50 repls
      wtreeqmc_wn_n2 : 7.08 (fn avg) 0.07 (rf avg) - 50 repls

400 bp, 1000 genes, bs support
           wastrid_s : 5.06 (fn avg) 0.05 (rf avg) - 50 repls
         aster_h_t16 : 5.04 (fn avg) 0.05 (rf avg) - 50 repls
      wtreeqmc_wh_n2 : 4.70 (fn avg) 0.05 (rf avg) - 50 repls
      wtreeqmc_wh_n1 : 4.96 (fn avg) 0.05 (rf avg) - 50 repls
      wtreeqmc_wh_n0 : 4.96 (fn avg) 0.05 (rf avg) - 50 repls
      wtreeqmc_ws_n2 : 4.96 (fn avg) 0.05 (rf avg) - 50 repls
      wtreeqmc_wn_n2 : 6.22 (fn avg) 0.06 (rf avg) - 50 repls

400 bp, 1000 genes, abayes support
           wastrid_s : 5.02 (fn avg) 0.05 (rf avg) - 50 repls
         aster_h_t16 : 5.06 (fn avg) 0.05 (rf avg) - 50 repls
      wtreeqmc_wh_n2 : 4.56 (fn avg) 0.05 (rf avg) - 50 repls
      wtreeqmc_wh_n1 : 5.14 (fn avg) 0.05 (rf avg) - 50 repls
      wtreeqmc_wh_n0 : 5.14 (fn avg) 0.05 (rf avg) - 50 repls
      wtreeqmc_ws_n2 : 5.14 (fn avg) 0.05 (rf avg) - 50 repls
      wtreeqmc_wn_n2 : 6.24 (fn avg) 0.06 (rf avg) - 50 repls

800 bp, 50 genes, bs support
           wastrid_s : 11.12 (fn avg) 0.11 (rf avg) - 50 repls
         aster_h_t16 : 10.18 (fn avg) 0.10 (rf avg) - 50 repls
      wtreeqmc_wh_n2 : 10.10 (fn avg) 0.10 (rf avg) - 50 repls
      wtreeqmc_wh_n1 : 10.32 (fn avg) 0.10 (rf avg) - 50 repls
      wtreeqmc_wh_n0 : 10.32 (fn avg) 0.10 (rf avg) - 50 repls
      wtreeqmc_ws_n2 : 10.32 (fn avg) 0.10 (rf avg) - 50 repls
      wtreeqmc_wn_n2 : 11.24 (fn avg) 0.12 (rf avg) - 50 repls

800 bp, 50 genes, abayes support
           wastrid_s : 10.92 (fn avg) 0.11 (rf avg) - 50 repls
         aster_h_t16 : 10.62 (fn avg) 0.11 (rf avg) - 50 repls
      wtreeqmc_wh_n2 : 9.48 (fn avg) 0.10 (rf avg) - 50 repls
      wtreeqmc_wh_n1 : 9.86 (fn avg) 0.10 (rf avg) - 50 repls
      wtreeqmc_wh_n0 : 9.86 (fn avg) 0.10 (rf avg) - 50 repls
      wtreeqmc_ws_n2 : 9.86 (fn avg) 0.10 (rf avg) - 50 repls
      wtreeqmc_wn_n2 : 11.26 (fn avg) 0.12 (rf avg) - 50 repls

800 bp, 200 genes, bs support
           wastrid_s : 6.76 (fn avg) 0.07 (rf avg) - 50 repls
         aster_h_t16 : 6.34 (fn avg) 0.06 (rf avg) - 50 repls
      wtreeqmc_wh_n2 : 6.16 (fn avg) 0.06 (rf avg) - 50 repls
      wtreeqmc_wh_n1 : 6.28 (fn avg) 0.06 (rf avg) - 50 repls
      wtreeqmc_wh_n0 : 6.28 (fn avg) 0.06 (rf avg) - 50 repls
      wtreeqmc_ws_n2 : 6.28 (fn avg) 0.06 (rf avg) - 50 repls
      wtreeqmc_wn_n2 : 6.86 (fn avg) 0.07 (rf avg) - 50 repls

800 bp, 200 genes, abayes support
           wastrid_s : 6.64 (fn avg) 0.07 (rf avg) - 50 repls
         aster_h_t16 : 6.28 (fn avg) 0.06 (rf avg) - 50 repls
      wtreeqmc_wh_n2 : 6.06 (fn avg) 0.06 (rf avg) - 50 repls
      wtreeqmc_wh_n1 : 6.12 (fn avg) 0.06 (rf avg) - 50 repls
      wtreeqmc_wh_n0 : 6.12 (fn avg) 0.06 (rf avg) - 50 repls
      wtreeqmc_ws_n2 : 6.12 (fn avg) 0.06 (rf avg) - 50 repls
      wtreeqmc_wn_n2 : 6.88 (fn avg) 0.07 (rf avg) - 50 repls

800 bp, 500 genes, bs support
           wastrid_s : 4.98 (fn avg) 0.05 (rf avg) - 50 repls
         aster_h_t16 : 4.80 (fn avg) 0.05 (rf avg) - 50 repls
      wtreeqmc_wh_n2 : 4.70 (fn avg) 0.05 (rf avg) - 50 repls
      wtreeqmc_wh_n1 : 4.82 (fn avg) 0.05 (rf avg) - 50 repls
      wtreeqmc_wh_n0 : 4.82 (fn avg) 0.05 (rf avg) - 50 repls
      wtreeqmc_ws_n2 : 4.82 (fn avg) 0.05 (rf avg) - 50 repls
      wtreeqmc_wn_n2 : 5.36 (fn avg) 0.06 (rf avg) - 50 repls

800 bp, 500 genes, abayes support
           wastrid_s : 4.78 (fn avg) 0.05 (rf avg) - 50 repls
         aster_h_t16 : 4.60 (fn avg) 0.05 (rf avg) - 50 repls
      wtreeqmc_wh_n2 : 4.28 (fn avg) 0.04 (rf avg) - 50 repls
      wtreeqmc_wh_n1 : 4.72 (fn avg) 0.05 (rf avg) - 50 repls
      wtreeqmc_wh_n0 : 4.72 (fn avg) 0.05 (rf avg) - 50 repls
      wtreeqmc_ws_n2 : 4.72 (fn avg) 0.05 (rf avg) - 50 repls
      wtreeqmc_wn_n2 : 5.34 (fn avg) 0.05 (rf avg) - 50 repls

800 bp, 1000 genes, bs support
           wastrid_s : 4.12 (fn avg) 0.04 (rf avg) - 50 repls
         aster_h_t16 : 4.10 (fn avg) 0.04 (rf avg) - 50 repls
      wtreeqmc_wh_n2 : 3.82 (fn avg) 0.04 (rf avg) - 50 repls
      wtreeqmc_wh_n1 : 4.04 (fn avg) 0.04 (rf avg) - 50 repls
      wtreeqmc_wh_n0 : 4.04 (fn avg) 0.04 (rf avg) - 50 repls
      wtreeqmc_ws_n2 : 4.04 (fn avg) 0.04 (rf avg) - 50 repls
      wtreeqmc_wn_n2 : 4.94 (fn avg) 0.05 (rf avg) - 50 repls

800 bp, 1000 genes, abayes support
           wastrid_s : 4.14 (fn avg) 0.04 (rf avg) - 50 repls
         aster_h_t16 : 4.10 (fn avg) 0.04 (rf avg) - 50 repls
      wtreeqmc_wh_n2 : 3.56 (fn avg) 0.04 (rf avg) - 50 repls
      wtreeqmc_wh_n1 : 4.06 (fn avg) 0.04 (rf avg) - 50 repls
      wtreeqmc_wh_n0 : 4.06 (fn avg) 0.04 (rf avg) - 50 repls
      wtreeqmc_ws_n2 : 4.06 (fn avg) 0.04 (rf avg) - 50 repls
      wtreeqmc_wn_n2 : 4.94 (fn avg) 0.05 (rf avg) - 50 repls

1600 bp, 50 genes, bs support
           wastrid_s : 10.14 (fn avg) 0.10 (rf avg) - 50 repls
         aster_h_t16 : 9.40 (fn avg) 0.10 (rf avg) - 50 repls
      wtreeqmc_wh_n2 : 9.10 (fn avg) 0.09 (rf avg) - 50 repls
      wtreeqmc_wh_n1 : 9.50 (fn avg) 0.10 (rf avg) - 50 repls
      wtreeqmc_wh_n0 : 9.50 (fn avg) 0.10 (rf avg) - 50 repls
      wtreeqmc_ws_n2 : 9.50 (fn avg) 0.10 (rf avg) - 50 repls
      wtreeqmc_wn_n2 : 9.52 (fn avg) 0.10 (rf avg) - 50 repls

1600 bp, 50 genes, abayes support
           wastrid_s : 10.24 (fn avg) 0.10 (rf avg) - 50 repls
         aster_h_t16 : 9.36 (fn avg) 0.10 (rf avg) - 50 repls
      wtreeqmc_wh_n2 : 8.46 (fn avg) 0.09 (rf avg) - 50 repls
      wtreeqmc_wh_n1 : 8.64 (fn avg) 0.09 (rf avg) - 50 repls
      wtreeqmc_wh_n0 : 8.64 (fn avg) 0.09 (rf avg) - 50 repls
      wtreeqmc_ws_n2 : 8.64 (fn avg) 0.09 (rf avg) - 50 repls
      wtreeqmc_wn_n2 : 9.52 (fn avg) 0.10 (rf avg) - 50 repls

1600 bp, 200 genes, bs support
           wastrid_s : 6.02 (fn avg) 0.06 (rf avg) - 50 repls
         aster_h_t16 : 5.70 (fn avg) 0.06 (rf avg) - 50 repls
      wtreeqmc_wh_n2 : 5.44 (fn avg) 0.06 (rf avg) - 50 repls
      wtreeqmc_wh_n1 : 5.50 (fn avg) 0.06 (rf avg) - 50 repls
      wtreeqmc_wh_n0 : 5.50 (fn avg) 0.06 (rf avg) - 50 repls
      wtreeqmc_ws_n2 : 5.50 (fn avg) 0.06 (rf avg) - 50 repls
      wtreeqmc_wn_n2 : 6.28 (fn avg) 0.06 (rf avg) - 50 repls

1600 bp, 200 genes, abayes support
           wastrid_s : 5.76 (fn avg) 0.06 (rf avg) - 50 repls
         aster_h_t16 : 5.54 (fn avg) 0.06 (rf avg) - 50 repls
      wtreeqmc_wh_n2 : 5.06 (fn avg) 0.05 (rf avg) - 50 repls
      wtreeqmc_wh_n1 : 5.26 (fn avg) 0.05 (rf avg) - 50 repls
      wtreeqmc_wh_n0 : 5.26 (fn avg) 0.05 (rf avg) - 50 repls
      wtreeqmc_ws_n2 : 5.26 (fn avg) 0.05 (rf avg) - 50 repls
      wtreeqmc_wn_n2 : 6.28 (fn avg) 0.06 (rf avg) - 50 repls

1600 bp, 500 genes, bs support
           wastrid_s : 4.32 (fn avg) 0.04 (rf avg) - 50 repls
         aster_h_t16 : 4.20 (fn avg) 0.04 (rf avg) - 50 repls
      wtreeqmc_wh_n2 : 3.98 (fn avg) 0.04 (rf avg) - 50 repls
      wtreeqmc_wh_n1 : 4.16 (fn avg) 0.04 (rf avg) - 50 repls
      wtreeqmc_wh_n0 : 4.16 (fn avg) 0.04 (rf avg) - 50 repls
      wtreeqmc_ws_n2 : 4.16 (fn avg) 0.04 (rf avg) - 50 repls
      wtreeqmc_wn_n2 : 4.80 (fn avg) 0.05 (rf avg) - 50 repls

1600 bp, 500 genes, abayes support
           wastrid_s : 4.00 (fn avg) 0.04 (rf avg) - 50 repls
         aster_h_t16 : 3.84 (fn avg) 0.04 (rf avg) - 50 repls
      wtreeqmc_wh_n2 : 3.46 (fn avg) 0.04 (rf avg) - 50 repls
      wtreeqmc_wh_n1 : 3.74 (fn avg) 0.04 (rf avg) - 50 repls
      wtreeqmc_wh_n0 : 3.74 (fn avg) 0.04 (rf avg) - 50 repls
      wtreeqmc_ws_n2 : 3.74 (fn avg) 0.04 (rf avg) - 50 repls
      wtreeqmc_wn_n2 : 4.80 (fn avg) 0.05 (rf avg) - 50 repls

1600 bp, 1000 genes, bs support
           wastrid_s : 3.52 (fn avg) 0.04 (rf avg) - 50 repls
         aster_h_t16 : 3.54 (fn avg) 0.04 (rf avg) - 50 repls
      wtreeqmc_wh_n2 : 3.44 (fn avg) 0.04 (rf avg) - 50 repls
      wtreeqmc_wh_n1 : 3.46 (fn avg) 0.04 (rf avg) - 50 repls
      wtreeqmc_wh_n0 : 3.46 (fn avg) 0.04 (rf avg) - 50 repls
      wtreeqmc_ws_n2 : 3.46 (fn avg) 0.04 (rf avg) - 50 repls
      wtreeqmc_wn_n2 : 4.02 (fn avg) 0.04 (rf avg) - 50 repls

1600 bp, 1000 genes, abayes support
           wastrid_s : 3.66 (fn avg) 0.04 (rf avg) - 50 repls
         aster_h_t16 : 3.46 (fn avg) 0.04 (rf avg) - 50 repls
      wtreeqmc_wh_n2 : 2.82 (fn avg) 0.03 (rf avg) - 50 repls
      wtreeqmc_wh_n1 : 3.34 (fn avg) 0.03 (rf avg) - 50 repls
      wtreeqmc_wh_n0 : 3.34 (fn avg) 0.03 (rf avg) - 50 repls
      wtreeqmc_ws_n2 : 3.34 (fn avg) 0.03 (rf avg) - 50 repls
      wtreeqmc_wn_n2 : 4.02 (fn avg) 0.04 (rf avg) - 50 repls
"""
