import pandas
import numpy
import sys

### DATA FROM FIGURE 2A in WASTRAL paper
### ALSO SEE FIGURE 2 in WASTRID paper

### IMPORTANT: We ran Asteroid and found that it had 100% error!
###            We looked at a few cases and found that the error was much better when using
###            -n option, which turns off the missing data correction. However, this is the
###            same as running ASTRID, and wASTRID outperforms ASTRID in these data.
###            Therefore, we didn't run Asteroid and don't include the results.

df = pandas.read_csv("data-all-error.csv", na_values='NA', keep_default_na=False)

mthds = ["ASTRID-ws",
         "ASTER-wh",
         "TQMC-wh_n2",
         "TQMC-wh_n1",
         "TQMC-wh_n0",
         "TQMC-ws_n2",
         "TQMC-wn_n2",
         "TQMC-n2"]

names = [" ASTRID-ws",
         "  ASTER-wh",
         "TQMC-wh_n2",
         "TQMC-wh_n1",
         "TQMC-wh_n0",
         "TQMC-ws_n2",
         "TQMC-wn_n2",
         "   TQMC-n2"]

nbpss = [200, 400, 800, 1600]
ngens = [50, 200, 500, 1000]
supps = ["bs", "abayes"]


for nbps in nbpss:
    for ngen in ngens:
        sys.stdout.write("%d bp, %d genes\n" % (nbps, ngen))
        for name, mthd in zip(names, mthds):
            sys.stdout.write("  %s : " % name)
            xdf = df[(df["NBPS"] == nbps) & 
                     (df["NGEN"] == ngen) &
                     (df["MTHD"] == mthd)]
            for supp in supps:
                ydf = xdf[xdf["SUPP"] == supp]
                data = ydf.SERF.values
                if data.size != 50:
                    sys.stdout.write("ERROR - wrong number of replicates!\n")

                nrepl = numpy.sum(~numpy.isnan(data))

                fnavg = numpy.mean(ydf.SEFN.values)
                fnavg = numpy.round(numpy.round(fnavg, 3), 2)

                rfavg = numpy.mean(ydf.SERF.values)
                rfavg = numpy.round(numpy.round(rfavg, 5), 4)

                sys.stdout.write("%1.4f (%s) " % (rfavg, supp))

            sys.stdout.write("\n")
        sys.stdout.write("\n")

"""
200 bp, 50 genes
   ASTRID-ws : 0.1631 (bs) 0.1504 (abayes) 
    ASTER-wh : 0.1553 (bs) 0.1457 (abayes) 
  TQMC-wh_n2 : 0.1539 (bs) 0.1414 (abayes) 
  TQMC-wh_n1 : 0.1563 (bs) 0.1457 (abayes) 
  TQMC-wh_n0 : 0.1674 (bs) 0.1586 (abayes) 
  TQMC-ws_n2 : 0.1561 (bs) 0.1455 (abayes) 
  TQMC-wn_n2 : 0.1720 (bs) 0.1739 (abayes) 
     TQMC-n2 : 0.1759 (bs) 0.1739 (abayes) 

200 bp, 200 genes
   ASTRID-ws : 0.0971 (bs) 0.0947 (abayes) 
    ASTER-wh : 0.0963 (bs) 0.0922 (abayes) 
  TQMC-wh_n2 : 0.0949 (bs) 0.0831 (abayes) 
  TQMC-wh_n1 : 0.0988 (bs) 0.0871 (abayes) 
  TQMC-wh_n0 : 0.1002 (bs) 0.0980 (abayes) 
  TQMC-ws_n2 : 0.0982 (bs) 0.0943 (abayes) 
  TQMC-wn_n2 : 0.1094 (bs) 0.1114 (abayes) 
     TQMC-n2 : 0.1096 (bs) 0.1114 (abayes) 

200 bp, 500 genes
   ASTRID-ws : 0.0774 (bs) 0.0753 (abayes) 
    ASTER-wh : 0.0753 (bs) 0.0716 (abayes) 
  TQMC-wh_n2 : 0.0755 (bs) 0.0692 (abayes) 
  TQMC-wh_n1 : 0.0792 (bs) 0.0729 (abayes) 
  TQMC-wh_n0 : 0.0802 (bs) 0.0778 (abayes) 
  TQMC-ws_n2 : 0.0794 (bs) 0.0784 (abayes) 
  TQMC-wn_n2 : 0.0929 (bs) 0.0924 (abayes) 
     TQMC-n2 : 0.0926 (bs) 0.0924 (abayes) 

200 bp, 1000 genes
   ASTRID-ws : 0.0669 (bs) 0.0622 (abayes) 
    ASTER-wh : 0.0653 (bs) 0.0631 (abayes) 
  TQMC-wh_n2 : 0.0673 (bs) 0.0535 (abayes) 
  TQMC-wh_n1 : 0.0716 (bs) 0.0610 (abayes) 
  TQMC-wh_n0 : 0.0702 (bs) 0.0684 (abayes) 
  TQMC-ws_n2 : 0.0704 (bs) 0.0592 (abayes) 
  TQMC-wn_n2 : 0.0796 (bs) 0.0800 (abayes) 
     TQMC-n2 : 0.0800 (bs) 0.0800 (abayes) 

400 bp, 50 genes
   ASTRID-ws : 0.1312 (bs) 0.1239 (abayes) 
    ASTER-wh : 0.1255 (bs) 0.1218 (abayes) 
  TQMC-wh_n2 : 0.1222 (bs) 0.1112 (abayes) 
  TQMC-wh_n1 : 0.1235 (bs) 0.1169 (abayes) 
  TQMC-wh_n0 : 0.1320 (bs) 0.1233 (abayes) 
  TQMC-ws_n2 : 0.1235 (bs) 0.1153 (abayes) 
  TQMC-wn_n2 : 0.1306 (bs) 0.1312 (abayes) 
     TQMC-n2 : 0.1310 (bs) 0.1312 (abayes) 

400 bp, 200 genes
   ASTRID-ws : 0.0796 (bs) 0.0743 (abayes) 
    ASTER-wh : 0.0767 (bs) 0.0722 (abayes) 
  TQMC-wh_n2 : 0.0729 (bs) 0.0671 (abayes) 
  TQMC-wh_n1 : 0.0769 (bs) 0.0712 (abayes) 
  TQMC-wh_n0 : 0.0833 (bs) 0.0778 (abayes) 
  TQMC-ws_n2 : 0.0747 (bs) 0.0718 (abayes) 
  TQMC-wn_n2 : 0.0853 (bs) 0.0859 (abayes) 
     TQMC-n2 : 0.0855 (bs) 0.0859 (abayes) 

400 bp, 500 genes
   ASTRID-ws : 0.0616 (bs) 0.0608 (abayes) 
    ASTER-wh : 0.0596 (bs) 0.0576 (abayes) 
  TQMC-wh_n2 : 0.0555 (bs) 0.0586 (abayes) 
  TQMC-wh_n1 : 0.0596 (bs) 0.0584 (abayes) 
  TQMC-wh_n0 : 0.0645 (bs) 0.0608 (abayes) 
  TQMC-ws_n2 : 0.0600 (bs) 0.0616 (abayes) 
  TQMC-wn_n2 : 0.0718 (bs) 0.0722 (abayes) 
     TQMC-n2 : 0.0722 (bs) 0.0722 (abayes) 

400 bp, 1000 genes
   ASTRID-ws : 0.0516 (bs) 0.0512 (abayes) 
    ASTER-wh : 0.0514 (bs) 0.0516 (abayes) 
  TQMC-wh_n2 : 0.0480 (bs) 0.0465 (abayes) 
  TQMC-wh_n1 : 0.0502 (bs) 0.0516 (abayes) 
  TQMC-wh_n0 : 0.0535 (bs) 0.0561 (abayes) 
  TQMC-ws_n2 : 0.0506 (bs) 0.0524 (abayes) 
  TQMC-wn_n2 : 0.0635 (bs) 0.0637 (abayes) 
     TQMC-n2 : 0.0633 (bs) 0.0637 (abayes) 

800 bp, 50 genes
   ASTRID-ws : 0.1135 (bs) 0.1114 (abayes) 
    ASTER-wh : 0.1039 (bs) 0.1084 (abayes) 
  TQMC-wh_n2 : 0.1031 (bs) 0.0967 (abayes) 
  TQMC-wh_n1 : 0.1043 (bs) 0.1037 (abayes) 
  TQMC-wh_n0 : 0.1104 (bs) 0.1126 (abayes) 
  TQMC-ws_n2 : 0.1053 (bs) 0.1006 (abayes) 
  TQMC-wn_n2 : 0.1147 (bs) 0.1149 (abayes) 
     TQMC-n2 : 0.1145 (bs) 0.1149 (abayes) 

800 bp, 200 genes
   ASTRID-ws : 0.0690 (bs) 0.0678 (abayes) 
    ASTER-wh : 0.0647 (bs) 0.0641 (abayes) 
  TQMC-wh_n2 : 0.0629 (bs) 0.0618 (abayes) 
  TQMC-wh_n1 : 0.0643 (bs) 0.0641 (abayes) 
  TQMC-wh_n0 : 0.0680 (bs) 0.0676 (abayes) 
  TQMC-ws_n2 : 0.0641 (bs) 0.0624 (abayes) 
  TQMC-wn_n2 : 0.0700 (bs) 0.0702 (abayes) 
     TQMC-n2 : 0.0698 (bs) 0.0702 (abayes) 

800 bp, 500 genes
   ASTRID-ws : 0.0508 (bs) 0.0488 (abayes) 
    ASTER-wh : 0.0490 (bs) 0.0469 (abayes) 
  TQMC-wh_n2 : 0.0480 (bs) 0.0437 (abayes) 
  TQMC-wh_n1 : 0.0496 (bs) 0.0459 (abayes) 
  TQMC-wh_n0 : 0.0535 (bs) 0.0496 (abayes) 
  TQMC-ws_n2 : 0.0492 (bs) 0.0482 (abayes) 
  TQMC-wn_n2 : 0.0547 (bs) 0.0545 (abayes) 
     TQMC-n2 : 0.0547 (bs) 0.0545 (abayes) 

800 bp, 1000 genes
   ASTRID-ws : 0.0420 (bs) 0.0422 (abayes) 
    ASTER-wh : 0.0418 (bs) 0.0418 (abayes) 
  TQMC-wh_n2 : 0.0390 (bs) 0.0363 (abayes) 
  TQMC-wh_n1 : 0.0426 (bs) 0.0406 (abayes) 
  TQMC-wh_n0 : 0.0449 (bs) 0.0465 (abayes) 
  TQMC-ws_n2 : 0.0412 (bs) 0.0414 (abayes) 
  TQMC-wn_n2 : 0.0504 (bs) 0.0504 (abayes) 
     TQMC-n2 : 0.0504 (bs) 0.0504 (abayes) 

1600 bp, 50 genes
   ASTRID-ws : 0.1035 (bs) 0.1045 (abayes) 
    ASTER-wh : 0.0959 (bs) 0.0955 (abayes) 
  TQMC-wh_n2 : 0.0929 (bs) 0.0863 (abayes) 
  TQMC-wh_n1 : 0.0953 (bs) 0.0931 (abayes) 
  TQMC-wh_n0 : 0.0992 (bs) 0.0994 (abayes) 
  TQMC-ws_n2 : 0.0969 (bs) 0.0882 (abayes) 
  TQMC-wn_n2 : 0.0971 (bs) 0.0971 (abayes) 
     TQMC-n2 : 0.0971 (bs) 0.0971 (abayes) 

1600 bp, 200 genes
   ASTRID-ws : 0.0614 (bs) 0.0588 (abayes) 
    ASTER-wh : 0.0582 (bs) 0.0565 (abayes) 
  TQMC-wh_n2 : 0.0555 (bs) 0.0516 (abayes) 
  TQMC-wh_n1 : 0.0588 (bs) 0.0549 (abayes) 
  TQMC-wh_n0 : 0.0616 (bs) 0.0580 (abayes) 
  TQMC-ws_n2 : 0.0561 (bs) 0.0537 (abayes) 
  TQMC-wn_n2 : 0.0641 (bs) 0.0641 (abayes) 
     TQMC-n2 : 0.0641 (bs) 0.0641 (abayes) 

1600 bp, 500 genes
   ASTRID-ws : 0.0441 (bs) 0.0408 (abayes) 
    ASTER-wh : 0.0429 (bs) 0.0392 (abayes) 
  TQMC-wh_n2 : 0.0406 (bs) 0.0353 (abayes) 
  TQMC-wh_n1 : 0.0426 (bs) 0.0386 (abayes) 
  TQMC-wh_n0 : 0.0443 (bs) 0.0422 (abayes) 
  TQMC-ws_n2 : 0.0424 (bs) 0.0382 (abayes) 
  TQMC-wn_n2 : 0.0490 (bs) 0.0490 (abayes) 
     TQMC-n2 : 0.0488 (bs) 0.0490 (abayes) 

1600 bp, 1000 genes
   ASTRID-ws : 0.0359 (bs) 0.0374 (abayes) 
    ASTER-wh : 0.0361 (bs) 0.0353 (abayes) 
  TQMC-wh_n2 : 0.0351 (bs) 0.0288 (abayes) 
  TQMC-wh_n1 : 0.0363 (bs) 0.0349 (abayes) 
  TQMC-wh_n0 : 0.0382 (bs) 0.0400 (abayes) 
  TQMC-ws_n2 : 0.0353 (bs) 0.0341 (abayes) 
  TQMC-wn_n2 : 0.0410 (bs) 0.0410 (abayes) 
     TQMC-n2 : 0.0410 (bs) 0.0410 (abayes)
"""
