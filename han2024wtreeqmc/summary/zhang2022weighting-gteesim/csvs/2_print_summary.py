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

                fnavg = numpy.mean(ydf.SEFN.values)
                fnavg = numpy.round(numpy.round(fnavg, 3), 2)

                rfavg = numpy.mean(ydf.SERF.values)
                rfavg = numpy.round(numpy.round(rfavg, 5), 4)

                print("    %s : %1.2f (fn avg) %1.4f (rf avg) - %d repls" % (name, fnavg, rfavg, nrepl))

            print()

"""
200 bp, 50 genes, bs support
           wastrid_s : 15.98 (fn avg) 0.1631 (rf avg) - 50 repls
         aster_h_t16 : 15.22 (fn avg) 0.1553 (rf avg) - 50 repls
      wtreeqmc_wh_n2 : 15.08 (fn avg) 0.1539 (rf avg) - 50 repls
      wtreeqmc_wh_n1 : 15.30 (fn avg) 0.1561 (rf avg) - 50 repls
      wtreeqmc_wh_n0 : 15.30 (fn avg) 0.1561 (rf avg) - 50 repls
      wtreeqmc_ws_n2 : 15.30 (fn avg) 0.1561 (rf avg) - 50 repls
      wtreeqmc_wn_n2 : 16.86 (fn avg) 0.1720 (rf avg) - 50 repls

200 bp, 50 genes, abayes support
           wastrid_s : 14.74 (fn avg) 0.1504 (rf avg) - 50 repls
         aster_h_t16 : 14.28 (fn avg) 0.1457 (rf avg) - 50 repls
      wtreeqmc_wh_n2 : 13.86 (fn avg) 0.1414 (rf avg) - 50 repls
      wtreeqmc_wh_n1 : 14.26 (fn avg) 0.1455 (rf avg) - 50 repls
      wtreeqmc_wh_n0 : 14.26 (fn avg) 0.1455 (rf avg) - 50 repls
      wtreeqmc_ws_n2 : 14.26 (fn avg) 0.1455 (rf avg) - 50 repls
      wtreeqmc_wn_n2 : 17.04 (fn avg) 0.1739 (rf avg) - 50 repls

200 bp, 200 genes, bs support
           wastrid_s : 9.52 (fn avg) 0.0971 (rf avg) - 50 repls
         aster_h_t16 : 9.44 (fn avg) 0.0963 (rf avg) - 50 repls
      wtreeqmc_wh_n2 : 9.30 (fn avg) 0.0949 (rf avg) - 50 repls
      wtreeqmc_wh_n1 : 9.62 (fn avg) 0.0982 (rf avg) - 50 repls
      wtreeqmc_wh_n0 : 9.62 (fn avg) 0.0982 (rf avg) - 50 repls
      wtreeqmc_ws_n2 : 9.62 (fn avg) 0.0982 (rf avg) - 50 repls
      wtreeqmc_wn_n2 : 10.72 (fn avg) 0.1094 (rf avg) - 50 repls

200 bp, 200 genes, abayes support
           wastrid_s : 9.28 (fn avg) 0.0947 (rf avg) - 50 repls
         aster_h_t16 : 9.04 (fn avg) 0.0922 (rf avg) - 50 repls
      wtreeqmc_wh_n2 : 8.14 (fn avg) 0.0831 (rf avg) - 50 repls
      wtreeqmc_wh_n1 : 9.24 (fn avg) 0.0943 (rf avg) - 50 repls
      wtreeqmc_wh_n0 : 9.24 (fn avg) 0.0943 (rf avg) - 50 repls
      wtreeqmc_ws_n2 : 9.24 (fn avg) 0.0943 (rf avg) - 50 repls
      wtreeqmc_wn_n2 : 10.92 (fn avg) 0.1114 (rf avg) - 50 repls

200 bp, 500 genes, bs support
           wastrid_s : 7.58 (fn avg) 0.0774 (rf avg) - 50 repls
         aster_h_t16 : 7.38 (fn avg) 0.0753 (rf avg) - 50 repls
      wtreeqmc_wh_n2 : 7.40 (fn avg) 0.0755 (rf avg) - 50 repls
      wtreeqmc_wh_n1 : 7.78 (fn avg) 0.0794 (rf avg) - 50 repls
      wtreeqmc_wh_n0 : 7.78 (fn avg) 0.0794 (rf avg) - 50 repls
      wtreeqmc_ws_n2 : 7.78 (fn avg) 0.0794 (rf avg) - 50 repls
      wtreeqmc_wn_n2 : 9.10 (fn avg) 0.0929 (rf avg) - 50 repls

200 bp, 500 genes, abayes support
           wastrid_s : 7.38 (fn avg) 0.0753 (rf avg) - 50 repls
         aster_h_t16 : 7.02 (fn avg) 0.0716 (rf avg) - 50 repls
      wtreeqmc_wh_n2 : 6.78 (fn avg) 0.0692 (rf avg) - 50 repls
      wtreeqmc_wh_n1 : 7.68 (fn avg) 0.0784 (rf avg) - 50 repls
      wtreeqmc_wh_n0 : 7.68 (fn avg) 0.0784 (rf avg) - 50 repls
      wtreeqmc_ws_n2 : 7.68 (fn avg) 0.0784 (rf avg) - 50 repls
      wtreeqmc_wn_n2 : 9.06 (fn avg) 0.0924 (rf avg) - 50 repls

200 bp, 1000 genes, bs support
           wastrid_s : 6.56 (fn avg) 0.0669 (rf avg) - 50 repls
         aster_h_t16 : 6.40 (fn avg) 0.0653 (rf avg) - 50 repls
      wtreeqmc_wh_n2 : 6.60 (fn avg) 0.0673 (rf avg) - 50 repls
      wtreeqmc_wh_n1 : 6.90 (fn avg) 0.0704 (rf avg) - 50 repls
      wtreeqmc_wh_n0 : 6.90 (fn avg) 0.0704 (rf avg) - 50 repls
      wtreeqmc_ws_n2 : 6.90 (fn avg) 0.0704 (rf avg) - 50 repls
      wtreeqmc_wn_n2 : 7.80 (fn avg) 0.0796 (rf avg) - 50 repls

200 bp, 1000 genes, abayes support
           wastrid_s : 6.10 (fn avg) 0.0622 (rf avg) - 50 repls
         aster_h_t16 : 6.18 (fn avg) 0.0631 (rf avg) - 50 repls
      wtreeqmc_wh_n2 : 5.24 (fn avg) 0.0535 (rf avg) - 50 repls
      wtreeqmc_wh_n1 : 5.80 (fn avg) 0.0592 (rf avg) - 50 repls
      wtreeqmc_wh_n0 : 5.80 (fn avg) 0.0592 (rf avg) - 50 repls
      wtreeqmc_ws_n2 : 5.80 (fn avg) 0.0592 (rf avg) - 50 repls
      wtreeqmc_wn_n2 : 7.84 (fn avg) 0.0800 (rf avg) - 50 repls

400 bp, 50 genes, bs support
           wastrid_s : 12.86 (fn avg) 0.1312 (rf avg) - 50 repls
         aster_h_t16 : 12.30 (fn avg) 0.1255 (rf avg) - 50 repls
      wtreeqmc_wh_n2 : 11.98 (fn avg) 0.1222 (rf avg) - 50 repls
      wtreeqmc_wh_n1 : 12.10 (fn avg) 0.1235 (rf avg) - 50 repls
      wtreeqmc_wh_n0 : 12.10 (fn avg) 0.1235 (rf avg) - 50 repls
      wtreeqmc_ws_n2 : 12.10 (fn avg) 0.1235 (rf avg) - 50 repls
      wtreeqmc_wn_n2 : 12.80 (fn avg) 0.1306 (rf avg) - 50 repls

400 bp, 50 genes, abayes support
           wastrid_s : 12.14 (fn avg) 0.1239 (rf avg) - 50 repls
         aster_h_t16 : 11.94 (fn avg) 0.1218 (rf avg) - 50 repls
      wtreeqmc_wh_n2 : 10.90 (fn avg) 0.1112 (rf avg) - 50 repls
      wtreeqmc_wh_n1 : 11.30 (fn avg) 0.1153 (rf avg) - 50 repls
      wtreeqmc_wh_n0 : 11.30 (fn avg) 0.1153 (rf avg) - 50 repls
      wtreeqmc_ws_n2 : 11.30 (fn avg) 0.1153 (rf avg) - 50 repls
      wtreeqmc_wn_n2 : 12.86 (fn avg) 0.1312 (rf avg) - 50 repls

400 bp, 200 genes, bs support
           wastrid_s : 7.80 (fn avg) 0.0796 (rf avg) - 50 repls
         aster_h_t16 : 7.52 (fn avg) 0.0767 (rf avg) - 50 repls
      wtreeqmc_wh_n2 : 7.14 (fn avg) 0.0729 (rf avg) - 50 repls
      wtreeqmc_wh_n1 : 7.32 (fn avg) 0.0747 (rf avg) - 50 repls
      wtreeqmc_wh_n0 : 7.32 (fn avg) 0.0747 (rf avg) - 50 repls
      wtreeqmc_ws_n2 : 7.32 (fn avg) 0.0747 (rf avg) - 50 repls
      wtreeqmc_wn_n2 : 8.36 (fn avg) 0.0853 (rf avg) - 50 repls

400 bp, 200 genes, abayes support
           wastrid_s : 7.28 (fn avg) 0.0743 (rf avg) - 50 repls
         aster_h_t16 : 7.08 (fn avg) 0.0722 (rf avg) - 50 repls
      wtreeqmc_wh_n2 : 6.58 (fn avg) 0.0671 (rf avg) - 50 repls
      wtreeqmc_wh_n1 : 7.04 (fn avg) 0.0718 (rf avg) - 50 repls
      wtreeqmc_wh_n0 : 7.04 (fn avg) 0.0718 (rf avg) - 50 repls
      wtreeqmc_ws_n2 : 7.04 (fn avg) 0.0718 (rf avg) - 50 repls
      wtreeqmc_wn_n2 : 8.42 (fn avg) 0.0859 (rf avg) - 50 repls

400 bp, 500 genes, bs support
           wastrid_s : 6.04 (fn avg) 0.0616 (rf avg) - 50 repls
         aster_h_t16 : 5.84 (fn avg) 0.0596 (rf avg) - 50 repls
      wtreeqmc_wh_n2 : 5.44 (fn avg) 0.0555 (rf avg) - 50 repls
      wtreeqmc_wh_n1 : 5.88 (fn avg) 0.0600 (rf avg) - 50 repls
      wtreeqmc_wh_n0 : 5.88 (fn avg) 0.0600 (rf avg) - 50 repls
      wtreeqmc_ws_n2 : 5.88 (fn avg) 0.0600 (rf avg) - 50 repls
      wtreeqmc_wn_n2 : 7.04 (fn avg) 0.0718 (rf avg) - 50 repls

400 bp, 500 genes, abayes support
           wastrid_s : 5.96 (fn avg) 0.0608 (rf avg) - 50 repls
         aster_h_t16 : 5.64 (fn avg) 0.0576 (rf avg) - 50 repls
      wtreeqmc_wh_n2 : 5.74 (fn avg) 0.0586 (rf avg) - 50 repls
      wtreeqmc_wh_n1 : 6.04 (fn avg) 0.0616 (rf avg) - 50 repls
      wtreeqmc_wh_n0 : 6.04 (fn avg) 0.0616 (rf avg) - 50 repls
      wtreeqmc_ws_n2 : 6.04 (fn avg) 0.0616 (rf avg) - 50 repls
      wtreeqmc_wn_n2 : 7.08 (fn avg) 0.0722 (rf avg) - 50 repls

400 bp, 1000 genes, bs support
           wastrid_s : 5.06 (fn avg) 0.0516 (rf avg) - 50 repls
         aster_h_t16 : 5.04 (fn avg) 0.0514 (rf avg) - 50 repls
      wtreeqmc_wh_n2 : 4.70 (fn avg) 0.0480 (rf avg) - 50 repls
      wtreeqmc_wh_n1 : 4.96 (fn avg) 0.0506 (rf avg) - 50 repls
      wtreeqmc_wh_n0 : 4.96 (fn avg) 0.0506 (rf avg) - 50 repls
      wtreeqmc_ws_n2 : 4.96 (fn avg) 0.0506 (rf avg) - 50 repls
      wtreeqmc_wn_n2 : 6.22 (fn avg) 0.0635 (rf avg) - 50 repls

400 bp, 1000 genes, abayes support
           wastrid_s : 5.02 (fn avg) 0.0512 (rf avg) - 50 repls
         aster_h_t16 : 5.06 (fn avg) 0.0516 (rf avg) - 50 repls
      wtreeqmc_wh_n2 : 4.56 (fn avg) 0.0465 (rf avg) - 50 repls
      wtreeqmc_wh_n1 : 5.14 (fn avg) 0.0524 (rf avg) - 50 repls
      wtreeqmc_wh_n0 : 5.14 (fn avg) 0.0524 (rf avg) - 50 repls
      wtreeqmc_ws_n2 : 5.14 (fn avg) 0.0524 (rf avg) - 50 repls
      wtreeqmc_wn_n2 : 6.24 (fn avg) 0.0637 (rf avg) - 50 repls

800 bp, 50 genes, bs support
           wastrid_s : 11.12 (fn avg) 0.1135 (rf avg) - 50 repls
         aster_h_t16 : 10.18 (fn avg) 0.1039 (rf avg) - 50 repls
      wtreeqmc_wh_n2 : 10.10 (fn avg) 0.1031 (rf avg) - 50 repls
      wtreeqmc_wh_n1 : 10.32 (fn avg) 0.1053 (rf avg) - 50 repls
      wtreeqmc_wh_n0 : 10.32 (fn avg) 0.1053 (rf avg) - 50 repls
      wtreeqmc_ws_n2 : 10.32 (fn avg) 0.1053 (rf avg) - 50 repls
      wtreeqmc_wn_n2 : 11.24 (fn avg) 0.1147 (rf avg) - 50 repls

800 bp, 50 genes, abayes support
           wastrid_s : 10.92 (fn avg) 0.1114 (rf avg) - 50 repls
         aster_h_t16 : 10.62 (fn avg) 0.1084 (rf avg) - 50 repls
      wtreeqmc_wh_n2 : 9.48 (fn avg) 0.0967 (rf avg) - 50 repls
      wtreeqmc_wh_n1 : 9.86 (fn avg) 0.1006 (rf avg) - 50 repls
      wtreeqmc_wh_n0 : 9.86 (fn avg) 0.1006 (rf avg) - 50 repls
      wtreeqmc_ws_n2 : 9.86 (fn avg) 0.1006 (rf avg) - 50 repls
      wtreeqmc_wn_n2 : 11.26 (fn avg) 0.1149 (rf avg) - 50 repls

800 bp, 200 genes, bs support
           wastrid_s : 6.76 (fn avg) 0.0690 (rf avg) - 50 repls
         aster_h_t16 : 6.34 (fn avg) 0.0647 (rf avg) - 50 repls
      wtreeqmc_wh_n2 : 6.16 (fn avg) 0.0629 (rf avg) - 50 repls
      wtreeqmc_wh_n1 : 6.28 (fn avg) 0.0641 (rf avg) - 50 repls
      wtreeqmc_wh_n0 : 6.28 (fn avg) 0.0641 (rf avg) - 50 repls
      wtreeqmc_ws_n2 : 6.28 (fn avg) 0.0641 (rf avg) - 50 repls
      wtreeqmc_wn_n2 : 6.86 (fn avg) 0.0700 (rf avg) - 50 repls

800 bp, 200 genes, abayes support
           wastrid_s : 6.64 (fn avg) 0.0678 (rf avg) - 50 repls
         aster_h_t16 : 6.28 (fn avg) 0.0641 (rf avg) - 50 repls
      wtreeqmc_wh_n2 : 6.06 (fn avg) 0.0618 (rf avg) - 50 repls
      wtreeqmc_wh_n1 : 6.12 (fn avg) 0.0624 (rf avg) - 50 repls
      wtreeqmc_wh_n0 : 6.12 (fn avg) 0.0624 (rf avg) - 50 repls
      wtreeqmc_ws_n2 : 6.12 (fn avg) 0.0624 (rf avg) - 50 repls
      wtreeqmc_wn_n2 : 6.88 (fn avg) 0.0702 (rf avg) - 50 repls

800 bp, 500 genes, bs support
           wastrid_s : 4.98 (fn avg) 0.0508 (rf avg) - 50 repls
         aster_h_t16 : 4.80 (fn avg) 0.0490 (rf avg) - 50 repls
      wtreeqmc_wh_n2 : 4.70 (fn avg) 0.0480 (rf avg) - 50 repls
      wtreeqmc_wh_n1 : 4.82 (fn avg) 0.0492 (rf avg) - 50 repls
      wtreeqmc_wh_n0 : 4.82 (fn avg) 0.0492 (rf avg) - 50 repls
      wtreeqmc_ws_n2 : 4.82 (fn avg) 0.0492 (rf avg) - 50 repls
      wtreeqmc_wn_n2 : 5.36 (fn avg) 0.0547 (rf avg) - 50 repls

800 bp, 500 genes, abayes support
           wastrid_s : 4.78 (fn avg) 0.0488 (rf avg) - 50 repls
         aster_h_t16 : 4.60 (fn avg) 0.0469 (rf avg) - 50 repls
      wtreeqmc_wh_n2 : 4.28 (fn avg) 0.0437 (rf avg) - 50 repls
      wtreeqmc_wh_n1 : 4.72 (fn avg) 0.0482 (rf avg) - 50 repls
      wtreeqmc_wh_n0 : 4.72 (fn avg) 0.0482 (rf avg) - 50 repls
      wtreeqmc_ws_n2 : 4.72 (fn avg) 0.0482 (rf avg) - 50 repls
      wtreeqmc_wn_n2 : 5.34 (fn avg) 0.0545 (rf avg) - 50 repls

800 bp, 1000 genes, bs support
           wastrid_s : 4.12 (fn avg) 0.0420 (rf avg) - 50 repls
         aster_h_t16 : 4.10 (fn avg) 0.0418 (rf avg) - 50 repls
      wtreeqmc_wh_n2 : 3.82 (fn avg) 0.0390 (rf avg) - 50 repls
      wtreeqmc_wh_n1 : 4.04 (fn avg) 0.0412 (rf avg) - 50 repls
      wtreeqmc_wh_n0 : 4.04 (fn avg) 0.0412 (rf avg) - 50 repls
      wtreeqmc_ws_n2 : 4.04 (fn avg) 0.0412 (rf avg) - 50 repls
      wtreeqmc_wn_n2 : 4.94 (fn avg) 0.0504 (rf avg) - 50 repls

800 bp, 1000 genes, abayes support
           wastrid_s : 4.14 (fn avg) 0.0422 (rf avg) - 50 repls
         aster_h_t16 : 4.10 (fn avg) 0.0418 (rf avg) - 50 repls
      wtreeqmc_wh_n2 : 3.56 (fn avg) 0.0363 (rf avg) - 50 repls
      wtreeqmc_wh_n1 : 4.06 (fn avg) 0.0414 (rf avg) - 50 repls
      wtreeqmc_wh_n0 : 4.06 (fn avg) 0.0414 (rf avg) - 50 repls
      wtreeqmc_ws_n2 : 4.06 (fn avg) 0.0414 (rf avg) - 50 repls
      wtreeqmc_wn_n2 : 4.94 (fn avg) 0.0504 (rf avg) - 50 repls

1600 bp, 50 genes, bs support
           wastrid_s : 10.14 (fn avg) 0.1035 (rf avg) - 50 repls
         aster_h_t16 : 9.40 (fn avg) 0.0959 (rf avg) - 50 repls
      wtreeqmc_wh_n2 : 9.10 (fn avg) 0.0929 (rf avg) - 50 repls
      wtreeqmc_wh_n1 : 9.50 (fn avg) 0.0969 (rf avg) - 50 repls
      wtreeqmc_wh_n0 : 9.50 (fn avg) 0.0969 (rf avg) - 50 repls
      wtreeqmc_ws_n2 : 9.50 (fn avg) 0.0969 (rf avg) - 50 repls
      wtreeqmc_wn_n2 : 9.52 (fn avg) 0.0971 (rf avg) - 50 repls

1600 bp, 50 genes, abayes support
           wastrid_s : 10.24 (fn avg) 0.1045 (rf avg) - 50 repls
         aster_h_t16 : 9.36 (fn avg) 0.0955 (rf avg) - 50 repls
      wtreeqmc_wh_n2 : 8.46 (fn avg) 0.0863 (rf avg) - 50 repls
      wtreeqmc_wh_n1 : 8.64 (fn avg) 0.0882 (rf avg) - 50 repls
      wtreeqmc_wh_n0 : 8.64 (fn avg) 0.0882 (rf avg) - 50 repls
      wtreeqmc_ws_n2 : 8.64 (fn avg) 0.0882 (rf avg) - 50 repls
      wtreeqmc_wn_n2 : 9.52 (fn avg) 0.0971 (rf avg) - 50 repls

1600 bp, 200 genes, bs support
           wastrid_s : 6.02 (fn avg) 0.0614 (rf avg) - 50 repls
         aster_h_t16 : 5.70 (fn avg) 0.0582 (rf avg) - 50 repls
      wtreeqmc_wh_n2 : 5.44 (fn avg) 0.0555 (rf avg) - 50 repls
      wtreeqmc_wh_n1 : 5.50 (fn avg) 0.0561 (rf avg) - 50 repls
      wtreeqmc_wh_n0 : 5.50 (fn avg) 0.0561 (rf avg) - 50 repls
      wtreeqmc_ws_n2 : 5.50 (fn avg) 0.0561 (rf avg) - 50 repls
      wtreeqmc_wn_n2 : 6.28 (fn avg) 0.0641 (rf avg) - 50 repls

1600 bp, 200 genes, abayes support
           wastrid_s : 5.76 (fn avg) 0.0588 (rf avg) - 50 repls
         aster_h_t16 : 5.54 (fn avg) 0.0565 (rf avg) - 50 repls
      wtreeqmc_wh_n2 : 5.06 (fn avg) 0.0516 (rf avg) - 50 repls
      wtreeqmc_wh_n1 : 5.26 (fn avg) 0.0537 (rf avg) - 50 repls
      wtreeqmc_wh_n0 : 5.26 (fn avg) 0.0537 (rf avg) - 50 repls
      wtreeqmc_ws_n2 : 5.26 (fn avg) 0.0537 (rf avg) - 50 repls
      wtreeqmc_wn_n2 : 6.28 (fn avg) 0.0641 (rf avg) - 50 repls

1600 bp, 500 genes, bs support
           wastrid_s : 4.32 (fn avg) 0.0441 (rf avg) - 50 repls
         aster_h_t16 : 4.20 (fn avg) 0.0429 (rf avg) - 50 repls
      wtreeqmc_wh_n2 : 3.98 (fn avg) 0.0406 (rf avg) - 50 repls
      wtreeqmc_wh_n1 : 4.16 (fn avg) 0.0424 (rf avg) - 50 repls
      wtreeqmc_wh_n0 : 4.16 (fn avg) 0.0424 (rf avg) - 50 repls
      wtreeqmc_ws_n2 : 4.16 (fn avg) 0.0424 (rf avg) - 50 repls
      wtreeqmc_wn_n2 : 4.80 (fn avg) 0.0490 (rf avg) - 50 repls

1600 bp, 500 genes, abayes support
           wastrid_s : 4.00 (fn avg) 0.0408 (rf avg) - 50 repls
         aster_h_t16 : 3.84 (fn avg) 0.0392 (rf avg) - 50 repls
      wtreeqmc_wh_n2 : 3.46 (fn avg) 0.0353 (rf avg) - 50 repls
      wtreeqmc_wh_n1 : 3.74 (fn avg) 0.0382 (rf avg) - 50 repls
      wtreeqmc_wh_n0 : 3.74 (fn avg) 0.0382 (rf avg) - 50 repls
      wtreeqmc_ws_n2 : 3.74 (fn avg) 0.0382 (rf avg) - 50 repls
      wtreeqmc_wn_n2 : 4.80 (fn avg) 0.0490 (rf avg) - 50 repls

1600 bp, 1000 genes, bs support
           wastrid_s : 3.52 (fn avg) 0.0359 (rf avg) - 50 repls
         aster_h_t16 : 3.54 (fn avg) 0.0361 (rf avg) - 50 repls
      wtreeqmc_wh_n2 : 3.44 (fn avg) 0.0351 (rf avg) - 50 repls
      wtreeqmc_wh_n1 : 3.46 (fn avg) 0.0353 (rf avg) - 50 repls
      wtreeqmc_wh_n0 : 3.46 (fn avg) 0.0353 (rf avg) - 50 repls
      wtreeqmc_ws_n2 : 3.46 (fn avg) 0.0353 (rf avg) - 50 repls
      wtreeqmc_wn_n2 : 4.02 (fn avg) 0.0410 (rf avg) - 50 repls

1600 bp, 1000 genes, abayes support
           wastrid_s : 3.66 (fn avg) 0.0374 (rf avg) - 50 repls
         aster_h_t16 : 3.46 (fn avg) 0.0353 (rf avg) - 50 repls
      wtreeqmc_wh_n2 : 2.82 (fn avg) 0.0288 (rf avg) - 50 repls
      wtreeqmc_wh_n1 : 3.34 (fn avg) 0.0341 (rf avg) - 50 repls
      wtreeqmc_wh_n0 : 3.34 (fn avg) 0.0341 (rf avg) - 50 repls
      wtreeqmc_ws_n2 : 3.34 (fn avg) 0.0341 (rf avg) - 50 repls
      wtreeqmc_wn_n2 : 4.02 (fn avg) 0.0410 (rf avg) - 50 repls
"""
