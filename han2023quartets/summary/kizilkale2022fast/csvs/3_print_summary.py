import pandas
import numpy
import sys

df = pandas.read_csv("data-all-error-and-timings_exfig4-6.csv", keep_default_na=False)

mthds = ["huntress_v0.1.2.0_default",
         "fastral",
         "treeqmcbip_v1.0.0_n2",
         "scistree_v1.2.0.6",
         "fastme_v2.1.5"]
names = ["  huntress", 
         "   fastral",
         "treeqmcbip",
         "  scistree",
         "    fastme"]

ncxnms = [["n1000", "m300"], ["n300", "m300"], ["n300", "m1000"]]
metrics = ["FN", "FP", "TIME", "NQS"]

for metric in metrics:
    print("Looking at %s" % metric)
    for ncxnm in ncxnms:
        [ncell, nmut] = ncxnm

        nrepl = 10

        print("Model : %s x %s" % (ncell, nmut))
        xdf = df[(df["NCELL"] == ncell) & (df["NMUT"] == nmut)]

        for mthd, name in zip(mthds, names):
            ydf = xdf[(xdf["MTHD"] == mthd)]

            if metric == "FN":
                data = ydf.SE_FN.values / ydf.NINT_TRUE.values
            elif metric == "FP":
                data = ydf.SE_FP.values / ydf.NINT_ESTI.values
            elif metric == "TIME":
                data = ydf.SECS.values / 60
            elif metric == "NQS":
                data = ydf.NQS.values

            if data.size != nrepl:
                sys.exit("ERROR - wrong number of replicates!")

            avg = numpy.mean(data)
            avg = numpy.round(numpy.round(avg, 5), 4)

            std = numpy.std(data)
            std = numpy.round(numpy.round(std, 5), 4)

            med = numpy.median(data)
            med = numpy.round(numpy.round(med, 5), 4)

            print("    %s : $%1.4f \pm %1.4f$ (avg) %1.4f (med)" % (name, avg, std, med))

    print()


"""
Looking at FN
Model : n1000 x m300
      huntress : $0.3845 \pm 0.0345$ (avg) 0.3788 (med)
       fastral : $0.6681 \pm 0.0416$ (avg) 0.6599 (med)
    treeqmcbip : $0.6224 \pm 0.0405$ (avg) 0.6111 (med)
      scistree : $0.6894 \pm 0.0559$ (avg) 0.6701 (med)
        fastme : $0.6408 \pm 0.0323$ (avg) 0.6414 (med)
Model : n300 x m300
      huntress : $0.2527 \pm 0.0422$ (avg) 0.2411 (med)
       fastral : $0.4648 \pm 0.0375$ (avg) 0.4620 (med)
    treeqmcbip : $0.3364 \pm 0.0480$ (avg) 0.3429 (med)
      scistree : $0.3731 \pm 0.0545$ (avg) 0.3679 (med)
        fastme : $0.3547 \pm 0.0486$ (avg) 0.3523 (med)
Model : n300 x m1000
      huntress : $0.0788 \pm 0.0159$ (avg) 0.0787 (med)
       fastral : $0.1092 \pm 0.0251$ (avg) 0.1053 (med)
    treeqmcbip : $0.0511 \pm 0.0226$ (avg) 0.0525 (med)
      scistree : $0.0313 \pm 0.0238$ (avg) 0.0235 (med)
        fastme : $0.0371 \pm 0.0124$ (avg) 0.0357 (med)

Looking at FP
Model : n1000 x m300
      huntress : $0.6993 \pm 0.0164$ (avg) 0.6966 (med)
       fastral : $0.8334 \pm 0.0203$ (avg) 0.8291 (med)
    treeqmcbip : $0.8089 \pm 0.0199$ (avg) 0.7999 (med)
      scistree : $0.8479 \pm 0.0278$ (avg) 0.8414 (med)
        fastme : $0.8288 \pm 0.0162$ (avg) 0.8290 (med)
Model : n300 x m300
      huntress : $0.5009 \pm 0.0371$ (avg) 0.5057 (med)
       fastral : $0.6495 \pm 0.0353$ (avg) 0.6530 (med)
    treeqmcbip : $0.5690 \pm 0.0310$ (avg) 0.5676 (med)
      scistree : $0.6085 \pm 0.0391$ (avg) 0.6056 (med)
        fastme : $0.5945 \pm 0.0374$ (avg) 0.5996 (med)
Model : n300 x m1000
      huntress : $0.5334 \pm 0.0132$ (avg) 0.5300 (med)
       fastral : $0.5536 \pm 0.0236$ (avg) 0.5493 (med)
    treeqmcbip : $0.5215 \pm 0.0320$ (avg) 0.5355 (med)
      scistree : $0.5649 \pm 0.0218$ (avg) 0.5659 (med)
        fastme : $0.5322 \pm 0.0221$ (avg) 0.5378 (med)

Looking at TIME
Model : n1000 x m300
      huntress : $3.1651 \pm 0.1335$ (avg) 3.2026 (med)
       fastral : $188.9739 \pm 38.9525$ (avg) 175.9366 (med)
    treeqmcbip : $114.5946 \pm 1.9964$ (avg) 113.7447 (med)
      scistree : $705.0790 \pm 87.2307$ (avg) 689.4959 (med)
        fastme : $2.5114 \pm 0.2255$ (avg) 2.5532 (med)
Model : n300 x m300
      huntress : $2.0495 \pm 0.1021$ (avg) 2.0013 (med)
       fastral : $4.8382 \pm 0.9492$ (avg) 4.5627 (med)
    treeqmcbip : $3.2482 \pm 0.0859$ (avg) 3.2056 (med)
      scistree : $16.3661 \pm 3.6336$ (avg) 16.6704 (med)
        fastme : $0.0402 \pm 0.0035$ (avg) 0.0399 (med)
Model : n300 x m1000
      huntress : $24.7036 \pm 1.8821$ (avg) 24.5651 (med)
       fastral : $7.6201 \pm 1.8643$ (avg) 6.9573 (med)
    treeqmcbip : $10.9433 \pm 0.3047$ (avg) 11.0146 (med)
      scistree : $53.4308 \pm 10.8458$ (avg) 51.0087 (med)
        fastme : $0.0371 \pm 0.0038$ (avg) 0.0376 (med)

Looking at NQS
Model : n1000 x m300
      huntress : $0.8622 \pm 0.0392$ (avg) 0.8602 (med)
       fastral : $0.8766 \pm 0.0340$ (avg) 0.8764 (med)
    treeqmcbip : $0.8802 \pm 0.0347$ (avg) 0.8807 (med)
      scistree : $0.8740 \pm 0.0376$ (avg) 0.8769 (med)
        fastme : $0.8766 \pm 0.0334$ (avg) 0.8781 (med)
Model : n300 x m300
      huntress : $0.8634 \pm 0.0402$ (avg) 0.8632 (med)
       fastral : $0.8815 \pm 0.0362$ (avg) 0.8838 (med)
    treeqmcbip : $0.8823 \pm 0.0354$ (avg) 0.8832 (med)
      scistree : $0.8811 \pm 0.0365$ (avg) 0.8841 (med)
        fastme : $0.8800 \pm 0.0368$ (avg) 0.8843 (med)
Model : n300 x m1000
      huntress : $0.8441 \pm 0.0853$ (avg) 0.8856 (med)
       fastral : $0.8606 \pm 0.0772$ (avg) 0.8958 (med)
    treeqmcbip : $0.8602 \pm 0.0773$ (avg) 0.8953 (med)
      scistree : $0.8598 \pm 0.0772$ (avg) 0.8952 (med)
        fastme : $0.8594 \pm 0.0770$ (avg) 0.8949 (med)
"""