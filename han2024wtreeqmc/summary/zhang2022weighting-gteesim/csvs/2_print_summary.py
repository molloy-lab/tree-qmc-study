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

df = pandas.read_csv("data-all-error-and-qscore.csv", na_values='NA', keep_default_na=False)

namemap = {}
namemap["ASTEROID"]   = "  ASTEROID"
namemap["ASTRID-ws"]  = " ASTRID-ws"
namemap["ASTER-wh"]   = "  ASTER-wh"
namemap["TQMC-wh_n2"] = "TQMC-wh_n2"
namemap["TQMC-wh_n1"] = "TQMC-wh_n1"
namemap["TQMC-wh_n0"] = "TQMC-wh_n0"
namemap["TQMC-ws_n2"] = "TQMC-ws_n2"
namemap["TQMC-wn_n2"] = "TQMC-wn_n2"
namemap["TQMC-n2"]    = "   TQMC-n2"

mthds = ["ASTEROID",
         "ASTRID-ws",
         "ASTER-wh",
         "TQMC-wh_n2",
         "TQMC-wh_n1",
         "TQMC-wh_n0",
         "TQMC-ws_n2",
         "TQMC-wn_n2",
         "TQMC-n2"]

nbpss = [200, 400, 800, 1600]
ngens = [50, 200, 500, 1000]



for nbps in nbpss:
    for ngen in ngens:
        sys.stdout.write("%d bp, %d genes\n" % (nbps, ngen))
        for mthd in mthds:
            if mthd == "ASTEROID":
               supps = ["bs"]
            else:
               supps = ["abayes", "bs"]

            sys.stdout.write("  %s : " % namemap[mthd])
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

                if (mthd == "TQMC-wh_n2") or (mthd == "ASTER-wh"):
                    qsavg = numpy.mean(ydf.QSCR.values)
                    qsavg = numpy.round(numpy.round(qsavg, 4), 3)

                    lppavg = numpy.nanmean(ydf.AVG_LPP.values)
                    lppavg = numpy.round(numpy.round(lppavg, 4), 3)

                    sys.stdout.write("    %1.4f (rf) %1.3f (qs) %1.3f (lpp) [%s]" % (rfavg, qsavg, lppavg, supp))
                else:
                    sys.stdout.write("    %1.4f (rf) [%s] " % (rfavg, supp))

            sys.stdout.write("\n")
        sys.stdout.write("\n")

"""
200 bp, 50 genes
    ASTEROID :     0.9971 (rf) [bs] 
   ASTRID-ws :     0.1504 (rf) [abayes]     0.1631 (rf) [bs] 
    ASTER-wh :     0.1457 (rf) 58551823.724 (qs) 0.882 (lpp) [abayes]    0.1553 (rf) 26134667.522 (qs) 0.895 (lpp) [bs]
  TQMC-wh_n2 :     0.1414 (rf) 58399244.157 (qs) 0.875 (lpp) [abayes]    0.1539 (rf) 26085601.984 (qs) 0.886 (lpp) [bs]
  TQMC-wh_n1 :     0.1457 (rf) [abayes]     0.1563 (rf) [bs] 
  TQMC-wh_n0 :     0.1586 (rf) [abayes]     0.1674 (rf) [bs] 
  TQMC-ws_n2 :     0.1455 (rf) [abayes]     0.1561 (rf) [bs] 
  TQMC-wn_n2 :     0.1739 (rf) [abayes]     0.1720 (rf) [bs] 
     TQMC-n2 :     0.1739 (rf) [abayes]     0.1759 (rf) [bs] 

200 bp, 200 genes
    ASTEROID :     0.9963 (rf) [bs] 
   ASTRID-ws :     0.0947 (rf) [abayes]     0.0971 (rf) [bs] 
    ASTER-wh :     0.0922 (rf) 233460448.566 (qs) 0.935 (lpp) [abayes]    0.0963 (rf) 103549926.193 (qs) 0.942 (lpp) [bs]
  TQMC-wh_n2 :     0.0831 (rf) 233292172.767 (qs) 0.928 (lpp) [abayes]    0.0949 (rf) 103472009.433 (qs) 0.935 (lpp) [bs]
  TQMC-wh_n1 :     0.0871 (rf) [abayes]     0.0988 (rf) [bs] 
  TQMC-wh_n0 :     0.0980 (rf) [abayes]     0.1002 (rf) [bs] 
  TQMC-ws_n2 :     0.0943 (rf) [abayes]     0.0982 (rf) [bs] 
  TQMC-wn_n2 :     0.1114 (rf) [abayes]     0.1094 (rf) [bs] 
     TQMC-n2 :     0.1114 (rf) [abayes]     0.1096 (rf) [bs] 

200 bp, 500 genes
    ASTEROID :     0.9969 (rf) [bs] 
   ASTRID-ws :     0.0753 (rf) [abayes]     0.0774 (rf) [bs] 
    ASTER-wh :     0.0716 (rf) 584176048.860 (qs) 0.957 (lpp) [abayes]    0.0753 (rf) 259997019.046 (qs) 0.963 (lpp) [bs]
  TQMC-wh_n2 :     0.0692 (rf) 583896058.814 (qs) 0.949 (lpp) [abayes]    0.0755 (rf) 259902733.394 (qs) 0.954 (lpp) [bs]
  TQMC-wh_n1 :     0.0729 (rf) [abayes]     0.0792 (rf) [bs] 
  TQMC-wh_n0 :     0.0778 (rf) [abayes]     0.0802 (rf) [bs] 
  TQMC-ws_n2 :     0.0784 (rf) [abayes]     0.0794 (rf) [bs] 
  TQMC-wn_n2 :     0.0924 (rf) [abayes]     0.0929 (rf) [bs] 
     TQMC-n2 :     0.0924 (rf) [abayes]     0.0926 (rf) [bs] 

200 bp, 1000 genes
    ASTEROID :     0.9967 (rf) [bs] 
   ASTRID-ws :     0.0622 (rf) [abayes]     0.0669 (rf) [bs] 
    ASTER-wh :     0.0631 (rf) 1166405478.746 (qs) 0.970 (lpp) [abayes]    0.0653 (rf) 519163476.734 (qs) 0.974 (lpp) [bs]
  TQMC-wh_n2 :     0.0535 (rf) 1166122490.996 (qs) 0.960 (lpp) [abayes]    0.0673 (rf) 519018402.664 (qs) 0.965 (lpp) [bs]
  TQMC-wh_n1 :     0.0610 (rf) [abayes]     0.0716 (rf) [bs] 
  TQMC-wh_n0 :     0.0684 (rf) [abayes]     0.0702 (rf) [bs] 
  TQMC-ws_n2 :     0.0592 (rf) [abayes]     0.0704 (rf) [bs] 
  TQMC-wn_n2 :     0.0800 (rf) [abayes]     0.0796 (rf) [bs] 
     TQMC-n2 :     0.0800 (rf) [abayes]     0.0800 (rf) [bs] 

400 bp, 50 genes
    ASTEROID :     0.9967 (rf) [bs] 
   ASTRID-ws :     0.1239 (rf) [abayes]     0.1312 (rf) [bs] 
    ASTER-wh :     0.1218 (rf) 69584349.708 (qs) 0.897 (lpp) [abayes]    0.1255 (rf) 41140105.400 (qs) 0.906 (lpp) [bs]
  TQMC-wh_n2 :     0.1112 (rf) 69490708.234 (qs) 0.892 (lpp) [abayes]    0.1222 (rf) 41093646.240 (qs) 0.900 (lpp) [bs]
  TQMC-wh_n1 :     0.1169 (rf) [abayes]     0.1235 (rf) [bs] 
  TQMC-wh_n0 :     0.1233 (rf) [abayes]     0.1320 (rf) [bs] 
  TQMC-ws_n2 :     0.1153 (rf) [abayes]     0.1235 (rf) [bs] 
  TQMC-wn_n2 :     0.1312 (rf) [abayes]     0.1306 (rf) [bs] 
     TQMC-n2 :     0.1312 (rf) [abayes]     0.1310 (rf) [bs] 

400 bp, 200 genes
    ASTEROID :     0.9967 (rf) [bs] 
   ASTRID-ws :     0.0743 (rf) [abayes]     0.0796 (rf) [bs] 
    ASTER-wh :     0.0722 (rf) 276678029.628 (qs) 0.947 (lpp) [abayes]    0.0767 (rf) 162556584.450 (qs) 0.949 (lpp) [bs]
  TQMC-wh_n2 :     0.0671 (rf) 276537166.200 (qs) 0.941 (lpp) [abayes]    0.0729 (rf) 162476539.587 (qs) 0.945 (lpp) [bs]
  TQMC-wh_n1 :     0.0712 (rf) [abayes]     0.0769 (rf) [bs] 
  TQMC-wh_n0 :     0.0778 (rf) [abayes]     0.0833 (rf) [bs] 
  TQMC-ws_n2 :     0.0718 (rf) [abayes]     0.0747 (rf) [bs] 
  TQMC-wn_n2 :     0.0859 (rf) [abayes]     0.0853 (rf) [bs] 
     TQMC-n2 :     0.0859 (rf) [abayes]     0.0855 (rf) [bs] 

400 bp, 500 genes
    ASTEROID :     0.9967 (rf) [bs] 
   ASTRID-ws :     0.0608 (rf) [abayes]     0.0616 (rf) [bs] 
    ASTER-wh :     0.0576 (rf) 693234057.498 (qs) 0.964 (lpp) [abayes]    0.0596 (rf) 408221017.368 (qs) 0.967 (lpp) [bs]
  TQMC-wh_n2 :     0.0586 (rf) 692889227.906 (qs) 0.958 (lpp) [abayes]    0.0555 (rf) 408114130.394 (qs) 0.962 (lpp) [bs]
  TQMC-wh_n1 :     0.0584 (rf) [abayes]     0.0596 (rf) [bs] 
  TQMC-wh_n0 :     0.0608 (rf) [abayes]     0.0645 (rf) [bs] 
  TQMC-ws_n2 :     0.0616 (rf) [abayes]     0.0600 (rf) [bs] 
  TQMC-wn_n2 :     0.0722 (rf) [abayes]     0.0718 (rf) [bs] 
     TQMC-n2 :     0.0722 (rf) [abayes]     0.0722 (rf) [bs] 

400 bp, 1000 genes
    ASTEROID :     0.9965 (rf) [bs] 
   ASTRID-ws :     0.0512 (rf) [abayes]     0.0516 (rf) [bs] 
    ASTER-wh :     0.0516 (rf) 1385134868.686 (qs) 0.975 (lpp) [abayes]    0.0514 (rf) 815341389.408 (qs) 0.976 (lpp) [bs]
  TQMC-wh_n2 :     0.0465 (rf) 1384713590.166 (qs) 0.968 (lpp) [abayes]    0.0480 (rf) 815204189.334 (qs) 0.972 (lpp) [bs]
  TQMC-wh_n1 :     0.0516 (rf) [abayes]     0.0502 (rf) [bs] 
  TQMC-wh_n0 :     0.0561 (rf) [abayes]     0.0535 (rf) [bs] 
  TQMC-ws_n2 :     0.0524 (rf) [abayes]     0.0506 (rf) [bs] 
  TQMC-wn_n2 :     0.0637 (rf) [abayes]     0.0635 (rf) [bs] 
     TQMC-n2 :     0.0637 (rf) [abayes]     0.0633 (rf) [bs] 

800 bp, 50 genes
    ASTEROID :     0.9969 (rf) [bs] 
   ASTRID-ws :     0.1114 (rf) [abayes]     0.1135 (rf) [bs] 
    ASTER-wh :     0.1084 (rf) 78432985.736 (qs) 0.907 (lpp) [abayes]    0.1039 (rf) 57402252.504 (qs) 0.910 (lpp) [bs]
  TQMC-wh_n2 :     0.0967 (rf) 78334173.614 (qs) 0.902 (lpp) [abayes]    0.1031 (rf) 57326298.984 (qs) 0.906 (lpp) [bs]
  TQMC-wh_n1 :     0.1037 (rf) [abayes]     0.1043 (rf) [bs] 
  TQMC-wh_n0 :     0.1126 (rf) [abayes]     0.1104 (rf) [bs] 
  TQMC-ws_n2 :     0.1006 (rf) [abayes]     0.1053 (rf) [bs] 
  TQMC-wn_n2 :     0.1149 (rf) [abayes]     0.1147 (rf) [bs] 
     TQMC-n2 :     0.1149 (rf) [abayes]     0.1145 (rf) [bs] 

800 bp, 200 genes
    ASTEROID :     0.9965 (rf) [bs] 
   ASTRID-ws :     0.0678 (rf) [abayes]     0.0690 (rf) [bs] 
    ASTER-wh :     0.0641 (rf) 311856335.696 (qs) 0.951 (lpp) [abayes]    0.0647 (rf) 226981945.849 (qs) 0.952 (lpp) [bs]
  TQMC-wh_n2 :     0.0618 (rf) 311704267.196 (qs) 0.946 (lpp) [abayes]    0.0629 (rf) 226886887.204 (qs) 0.948 (lpp) [bs]
  TQMC-wh_n1 :     0.0641 (rf) [abayes]     0.0643 (rf) [bs] 
  TQMC-wh_n0 :     0.0676 (rf) [abayes]     0.0680 (rf) [bs] 
  TQMC-ws_n2 :     0.0624 (rf) [abayes]     0.0641 (rf) [bs] 
  TQMC-wn_n2 :     0.0702 (rf) [abayes]     0.0700 (rf) [bs] 
     TQMC-n2 :     0.0702 (rf) [abayes]     0.0698 (rf) [bs] 

800 bp, 500 genes
    ASTEROID :     0.9967 (rf) [bs] 
   ASTRID-ws :     0.0488 (rf) [abayes]     0.0508 (rf) [bs] 
    ASTER-wh :     0.0469 (rf) 780980877.998 (qs) 0.968 (lpp) [abayes]    0.0490 (rf) 568984223.244 (qs) 0.969 (lpp) [bs]
  TQMC-wh_n2 :     0.0437 (rf) 780720362.602 (qs) 0.965 (lpp) [abayes]    0.0480 (rf) 568934647.792 (qs) 0.965 (lpp) [bs]
  TQMC-wh_n1 :     0.0459 (rf) [abayes]     0.0496 (rf) [bs] 
  TQMC-wh_n0 :     0.0496 (rf) [abayes]     0.0535 (rf) [bs] 
  TQMC-ws_n2 :     0.0482 (rf) [abayes]     0.0492 (rf) [bs] 
  TQMC-wn_n2 :     0.0545 (rf) [abayes]     0.0547 (rf) [bs] 
     TQMC-n2 :     0.0545 (rf) [abayes]     0.0547 (rf) [bs] 

800 bp, 1000 genes
    ASTEROID :     0.9965 (rf) [bs] 
   ASTRID-ws :     0.0422 (rf) [abayes]     0.0420 (rf) [bs] 
    ASTER-wh :     0.0418 (rf) 1561215385.824 (qs) 0.977 (lpp) [abayes]    0.0418 (rf) 1136796247.816 (qs) 0.977 (lpp) [bs]
  TQMC-wh_n2 :     0.0363 (rf) 1560763710.244 (qs) 0.972 (lpp) [abayes]    0.0390 (rf) 1136724808.972 (qs) 0.975 (lpp) [bs]
  TQMC-wh_n1 :     0.0406 (rf) [abayes]     0.0426 (rf) [bs] 
  TQMC-wh_n0 :     0.0465 (rf) [abayes]     0.0449 (rf) [bs] 
  TQMC-ws_n2 :     0.0414 (rf) [abayes]     0.0412 (rf) [bs] 
  TQMC-wn_n2 :     0.0504 (rf) [abayes]     0.0504 (rf) [bs] 
     TQMC-n2 :     0.0504 (rf) [abayes]     0.0504 (rf) [bs] 

1600 bp, 50 genes
    ASTEROID :     0.9965 (rf) [bs] 
   ASTRID-ws :     0.1045 (rf) [abayes]     0.1035 (rf) [bs] 
    ASTER-wh :     0.0955 (rf) 84367367.694 (qs) 0.910 (lpp) [abayes]    0.0959 (rf) 70544751.609 (qs) 0.913 (lpp) [bs]
  TQMC-wh_n2 :     0.0863 (rf) 84274386.534 (qs) 0.907 (lpp) [abayes]    0.0929 (rf) 70482876.296 (qs) 0.910 (lpp) [bs]
  TQMC-wh_n1 :     0.0931 (rf) [abayes]     0.0953 (rf) [bs] 
  TQMC-wh_n0 :     0.0994 (rf) [abayes]     0.0992 (rf) [bs] 
  TQMC-ws_n2 :     0.0882 (rf) [abayes]     0.0969 (rf) [bs] 
  TQMC-wn_n2 :     0.0971 (rf) [abayes]     0.0971 (rf) [bs] 
     TQMC-n2 :     0.0971 (rf) [abayes]     0.0971 (rf) [bs] 

1600 bp, 200 genes
    ASTEROID :     0.9963 (rf) [bs] 
   ASTRID-ws :     0.0588 (rf) [abayes]     0.0614 (rf) [bs] 
    ASTER-wh :     0.0565 (rf) 336838250.574 (qs) 0.953 (lpp) [abayes]    0.0582 (rf) 280404775.642 (qs) 0.954 (lpp) [bs]
  TQMC-wh_n2 :     0.0516 (rf) 336695383.764 (qs) 0.949 (lpp) [abayes]    0.0555 (rf) 280321251.486 (qs) 0.951 (lpp) [bs]
  TQMC-wh_n1 :     0.0549 (rf) [abayes]     0.0588 (rf) [bs] 
  TQMC-wh_n0 :     0.0580 (rf) [abayes]     0.0616 (rf) [bs] 
  TQMC-ws_n2 :     0.0537 (rf) [abayes]     0.0561 (rf) [bs] 
  TQMC-wn_n2 :     0.0641 (rf) [abayes]     0.0641 (rf) [bs] 
     TQMC-n2 :     0.0641 (rf) [abayes]     0.0641 (rf) [bs] 

1600 bp, 500 genes
    ASTEROID :     0.9967 (rf) [bs] 
   ASTRID-ws :     0.0408 (rf) [abayes]     0.0441 (rf) [bs] 
    ASTER-wh :     0.0392 (rf) 843689756.406 (qs) 0.968 (lpp) [abayes]    0.0429 (rf) 702530104.774 (qs) 0.970 (lpp) [bs]
  TQMC-wh_n2 :     0.0353 (rf) 843597247.792 (qs) 0.966 (lpp) [abayes]    0.0406 (rf) 702424823.536 (qs) 0.967 (lpp) [bs]
  TQMC-wh_n1 :     0.0386 (rf) [abayes]     0.0426 (rf) [bs] 
  TQMC-wh_n0 :     0.0422 (rf) [abayes]     0.0443 (rf) [bs] 
  TQMC-ws_n2 :     0.0382 (rf) [abayes]     0.0424 (rf) [bs] 
  TQMC-wn_n2 :     0.0490 (rf) [abayes]     0.0490 (rf) [bs] 
     TQMC-n2 :     0.0490 (rf) [abayes]     0.0488 (rf) [bs] 

1600 bp, 1000 genes
    ASTEROID :     0.9965 (rf) [bs] 
   ASTRID-ws :     0.0374 (rf) [abayes]     0.0359 (rf) [bs] 
    ASTER-wh :     0.0353 (rf) 1687283172.316 (qs) 0.978 (lpp) [abayes]    0.0361 (rf) 1404102591.908 (qs) 0.978 (lpp) [bs]
  TQMC-wh_n2 :     0.0288 (rf) 1687128469.576 (qs) 0.974 (lpp) [abayes]    0.0351 (rf) 1403984848.998 (qs) 0.975 (lpp) [bs]
  TQMC-wh_n1 :     0.0349 (rf) [abayes]     0.0363 (rf) [bs] 
  TQMC-wh_n0 :     0.0400 (rf) [abayes]     0.0382 (rf) [bs] 
  TQMC-ws_n2 :     0.0341 (rf) [abayes]     0.0353 (rf) [bs] 
  TQMC-wn_n2 :     0.0410 (rf) [abayes]     0.0410 (rf) [bs] 
     TQMC-n2 :     0.0410 (rf) [abayes]     0.0410 (rf) [bs] 
"""
