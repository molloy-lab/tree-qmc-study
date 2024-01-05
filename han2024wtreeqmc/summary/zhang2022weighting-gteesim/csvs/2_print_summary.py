import pandas
import numpy
import sys

#sys.exit("DONE RUNNING")

### DATA FROM FIGURE 2A in WASTRAL paper
### ALSO SEE FIGURE 2 in WASTRID paper

### IMPORTANT: We ran Asteroid and found that it had 100% error!
###            We looked at a few cases and found that the error was much better when using
###            -n option, which turns off the missing data correction. However, this is the
###            same as running ASTRID, and wASTRID outperforms ASTRID in these data.
###            Therefore, we didn't run Asteroid and don't include the results.

df = pandas.read_csv("data-all-error.csv", na_values='NA', keep_default_na=False)

mthds = ["wastrid_s",
         "aster_h_t16",
         "wtreeqmc_wh_n2",
         "wtreeqmc_wh_n1",
         "wtreeqmc_wh_n0",
         "wtreeqmc_ws_n2",
         "wtreeqmc_wn_n2"]

names = ["     wastrid_s",
         "   aster_h_t16",
         "wtreeqmc_wh_n2",
         "wtreeqmc_wh_n1",
         "wtreeqmc_wh_n0",
         "wtreeqmc_ws_n2",
         "wtreeqmc_wn_n2"]

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
                if mthd == "wtreeqmc_wn_n2":
                    if supp == "bs":
                        xsupp = "none_refinepoly"
                    else:
                        xsupp = "none_keeppoly"
                else:
                    xsupp = supp
                ydf = xdf[xdf["SUPP"] == xsupp]

                data = ydf.SERF.values
                #if data.size != 50:
                #    sys.stdout.write("ERROR - wrong number of replicates!\n")

                nrepl = numpy.sum(~numpy.isnan(data))

                fnavg = numpy.mean(ydf.SEFN.values)
                fnavg = numpy.round(numpy.round(fnavg, 3), 2)

                rfavg = numpy.mean(ydf.SERF.values)
                rfavg = numpy.round(numpy.round(rfavg, 5), 4)

                sys.stdout.write("%1.4f (%s) " % (rfavg, xsupp))

            sys.stdout.write("\n")
        sys.stdout.write("\n")

"""
200 bp, 50 genes
       wastrid_s : 0.1631 (bs) 0.1504 (abayes) 
     aster_h_t16 : 0.1553 (bs) 0.1457 (abayes) 
  wtreeqmc_wh_n2 : 0.1539 (bs) 0.1414 (abayes) 
  wtreeqmc_wh_n1 : 0.1561 (bs) 0.1455 (abayes) 
  wtreeqmc_wh_n0 : 0.1561 (bs) 0.1455 (abayes) 
  wtreeqmc_ws_n2 : 0.1561 (bs) 0.1455 (abayes) 
  wtreeqmc_wn_n2 : 0.1739 (none_refinepoly) 0.1720 (none_keeppoly) 

200 bp, 200 genes
       wastrid_s : 0.0971 (bs) 0.0947 (abayes) 
     aster_h_t16 : 0.0963 (bs) 0.0922 (abayes) 
  wtreeqmc_wh_n2 : 0.0949 (bs) 0.0831 (abayes) 
  wtreeqmc_wh_n1 : 0.0982 (bs) 0.0943 (abayes) 
  wtreeqmc_wh_n0 : 0.0982 (bs) 0.0943 (abayes) 
  wtreeqmc_ws_n2 : 0.0982 (bs) 0.0943 (abayes) 
  wtreeqmc_wn_n2 : 0.1114 (none_refinepoly) 0.1094 (none_keeppoly) 

200 bp, 500 genes
       wastrid_s : 0.0774 (bs) 0.0753 (abayes) 
     aster_h_t16 : 0.0753 (bs) 0.0716 (abayes) 
  wtreeqmc_wh_n2 : 0.0755 (bs) 0.0692 (abayes) 
  wtreeqmc_wh_n1 : 0.0794 (bs) 0.0784 (abayes) 
  wtreeqmc_wh_n0 : 0.0794 (bs) 0.0784 (abayes) 
  wtreeqmc_ws_n2 : 0.0794 (bs) 0.0784 (abayes) 
  wtreeqmc_wn_n2 : 0.0924 (none_refinepoly) 0.0929 (none_keeppoly) 

200 bp, 1000 genes
       wastrid_s : 0.0669 (bs) 0.0622 (abayes) 
     aster_h_t16 : 0.0653 (bs) 0.0631 (abayes) 
  wtreeqmc_wh_n2 : 0.0673 (bs) 0.0535 (abayes) 
  wtreeqmc_wh_n1 : 0.0704 (bs) 0.0592 (abayes) 
  wtreeqmc_wh_n0 : 0.0704 (bs) 0.0592 (abayes) 
  wtreeqmc_ws_n2 : 0.0704 (bs) 0.0592 (abayes) 
  wtreeqmc_wn_n2 : 0.0800 (none_refinepoly) 0.0796 (none_keeppoly) 

400 bp, 50 genes
       wastrid_s : 0.1312 (bs) 0.1239 (abayes) 
     aster_h_t16 : 0.1255 (bs) 0.1218 (abayes) 
  wtreeqmc_wh_n2 : 0.1222 (bs) 0.1112 (abayes) 
  wtreeqmc_wh_n1 : 0.1235 (bs) 0.1153 (abayes) 
  wtreeqmc_wh_n0 : 0.1235 (bs) 0.1153 (abayes) 
  wtreeqmc_ws_n2 : 0.1235 (bs) 0.1153 (abayes) 
  wtreeqmc_wn_n2 : 0.1312 (none_refinepoly) 0.1306 (none_keeppoly) 

400 bp, 200 genes
       wastrid_s : 0.0796 (bs) 0.0743 (abayes) 
     aster_h_t16 : 0.0767 (bs) 0.0722 (abayes) 
  wtreeqmc_wh_n2 : 0.0729 (bs) 0.0671 (abayes) 
  wtreeqmc_wh_n1 : 0.0747 (bs) 0.0718 (abayes) 
  wtreeqmc_wh_n0 : 0.0747 (bs) 0.0718 (abayes) 
  wtreeqmc_ws_n2 : 0.0747 (bs) 0.0718 (abayes) 
  wtreeqmc_wn_n2 : 0.0859 (none_refinepoly) 0.0853 (none_keeppoly) 

400 bp, 500 genes
       wastrid_s : 0.0616 (bs) 0.0608 (abayes) 
     aster_h_t16 : 0.0596 (bs) 0.0576 (abayes) 
  wtreeqmc_wh_n2 : 0.0555 (bs) 0.0586 (abayes) 
  wtreeqmc_wh_n1 : 0.0600 (bs) 0.0616 (abayes) 
  wtreeqmc_wh_n0 : 0.0600 (bs) 0.0616 (abayes) 
  wtreeqmc_ws_n2 : 0.0600 (bs) 0.0616 (abayes) 
  wtreeqmc_wn_n2 : 0.0722 (none_refinepoly) 0.0718 (none_keeppoly) 

400 bp, 1000 genes
       wastrid_s : 0.0516 (bs) 0.0512 (abayes) 
     aster_h_t16 : 0.0514 (bs) 0.0516 (abayes) 
  wtreeqmc_wh_n2 : 0.0480 (bs) 0.0465 (abayes) 
  wtreeqmc_wh_n1 : 0.0506 (bs) 0.0524 (abayes) 
  wtreeqmc_wh_n0 : 0.0506 (bs) 0.0524 (abayes) 
  wtreeqmc_ws_n2 : 0.0506 (bs) 0.0524 (abayes) 
  wtreeqmc_wn_n2 : 0.0637 (none_refinepoly) 0.0635 (none_keeppoly) 

800 bp, 50 genes
       wastrid_s : 0.1135 (bs) 0.1114 (abayes) 
     aster_h_t16 : 0.1039 (bs) 0.1084 (abayes) 
  wtreeqmc_wh_n2 : 0.1031 (bs) 0.0967 (abayes) 
  wtreeqmc_wh_n1 : 0.1053 (bs) 0.1006 (abayes) 
  wtreeqmc_wh_n0 : 0.1053 (bs) 0.1006 (abayes) 
  wtreeqmc_ws_n2 : 0.1053 (bs) 0.1006 (abayes) 
  wtreeqmc_wn_n2 : 0.1149 (none_refinepoly) 0.1147 (none_keeppoly) 

800 bp, 200 genes
       wastrid_s : 0.0690 (bs) 0.0678 (abayes) 
     aster_h_t16 : 0.0647 (bs) 0.0641 (abayes) 
  wtreeqmc_wh_n2 : 0.0629 (bs) 0.0618 (abayes) 
  wtreeqmc_wh_n1 : 0.0641 (bs) 0.0624 (abayes) 
  wtreeqmc_wh_n0 : 0.0641 (bs) 0.0624 (abayes) 
  wtreeqmc_ws_n2 : 0.0641 (bs) 0.0624 (abayes) 
  wtreeqmc_wn_n2 : 0.0702 (none_refinepoly) 0.0700 (none_keeppoly) 

800 bp, 500 genes
       wastrid_s : 0.0508 (bs) 0.0488 (abayes) 
     aster_h_t16 : 0.0490 (bs) 0.0469 (abayes) 
  wtreeqmc_wh_n2 : 0.0480 (bs) 0.0437 (abayes) 
  wtreeqmc_wh_n1 : 0.0492 (bs) 0.0482 (abayes) 
  wtreeqmc_wh_n0 : 0.0492 (bs) 0.0482 (abayes) 
  wtreeqmc_ws_n2 : 0.0492 (bs) 0.0482 (abayes) 
  wtreeqmc_wn_n2 : 0.0545 (none_refinepoly) 0.0547 (none_keeppoly) 

800 bp, 1000 genes
       wastrid_s : 0.0420 (bs) 0.0422 (abayes) 
     aster_h_t16 : 0.0418 (bs) 0.0418 (abayes) 
  wtreeqmc_wh_n2 : 0.0390 (bs) 0.0363 (abayes) 
  wtreeqmc_wh_n1 : 0.0412 (bs) 0.0414 (abayes) 
  wtreeqmc_wh_n0 : 0.0412 (bs) 0.0414 (abayes) 
  wtreeqmc_ws_n2 : 0.0412 (bs) 0.0414 (abayes) 
  wtreeqmc_wn_n2 : 0.0504 (none_refinepoly) 0.0504 (none_keeppoly) 

1600 bp, 50 genes
       wastrid_s : 0.1035 (bs) 0.1045 (abayes) 
     aster_h_t16 : 0.0959 (bs) 0.0955 (abayes) 
  wtreeqmc_wh_n2 : 0.0929 (bs) 0.0863 (abayes) 
  wtreeqmc_wh_n1 : 0.0969 (bs) 0.0882 (abayes) 
  wtreeqmc_wh_n0 : 0.0969 (bs) 0.0882 (abayes) 
  wtreeqmc_ws_n2 : 0.0969 (bs) 0.0882 (abayes) 
  wtreeqmc_wn_n2 : 0.0971 (none_refinepoly) 0.0971 (none_keeppoly) 

1600 bp, 200 genes
       wastrid_s : 0.0614 (bs) 0.0588 (abayes) 
     aster_h_t16 : 0.0582 (bs) 0.0565 (abayes) 
  wtreeqmc_wh_n2 : 0.0555 (bs) 0.0516 (abayes) 
  wtreeqmc_wh_n1 : 0.0561 (bs) 0.0537 (abayes) 
  wtreeqmc_wh_n0 : 0.0561 (bs) 0.0537 (abayes) 
  wtreeqmc_ws_n2 : 0.0561 (bs) 0.0537 (abayes) 
  wtreeqmc_wn_n2 : 0.0641 (none_refinepoly) 0.0641 (none_keeppoly) 

1600 bp, 500 genes
       wastrid_s : 0.0441 (bs) 0.0408 (abayes) 
     aster_h_t16 : 0.0429 (bs) 0.0392 (abayes) 
  wtreeqmc_wh_n2 : 0.0406 (bs) 0.0353 (abayes) 
  wtreeqmc_wh_n1 : 0.0424 (bs) 0.0382 (abayes) 
  wtreeqmc_wh_n0 : 0.0424 (bs) 0.0382 (abayes) 
  wtreeqmc_ws_n2 : 0.0424 (bs) 0.0382 (abayes) 
  wtreeqmc_wn_n2 : 0.0490 (none_refinepoly) 0.0490 (none_keeppoly) 

1600 bp, 1000 genes
       wastrid_s : 0.0359 (bs) 0.0374 (abayes) 
     aster_h_t16 : 0.0361 (bs) 0.0353 (abayes) 
  wtreeqmc_wh_n2 : 0.0351 (bs) 0.0288 (abayes) 
  wtreeqmc_wh_n1 : 0.0353 (bs) 0.0341 (abayes) 
  wtreeqmc_wh_n0 : 0.0353 (bs) 0.0341 (abayes) 
  wtreeqmc_ws_n2 : 0.0353 (bs) 0.0341 (abayes) 
  wtreeqmc_wn_n2 : 0.0410 (none_refinepoly) 0.0410 (none_keeppoly) 
"""
