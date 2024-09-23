import pandas
import numpy
import sys
import warnings

def report_stats(df, mthds, ngen, ntax, hght, rate):
    sys.stdout.write("%d & %s & %s & %d" % (ntax, hght, rate, ngen))

    nrepl = len(df[(df["MTHD"] == "CA-ML")].REPL.values)

    keep = []
    minval = 1
    for ind, mthd in enumerate(mthds):
        xdf = df[(df["MTHD"] == mthd)]
        vals = xdf.SEFNR.values

        if len(vals) != nrepl:
            sys.exit("\nERROR - inconsistent number of replicates")
        if len(vals) > 50:
            sys.exit("\nERROR - too many replicates")

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            fnravg = numpy.nanmean(vals)

        x = round(round(fnravg, 5), 4)
        keep.append(x)

        if x < minval:
            minval = x
            minind = []
            minind.append(ind)
        elif x == minval:
            minind.append(ind)

    minind = set(minind)
    for ind, x in enumerate(keep):
        if ind in minind:
            sys.stdout.write(" & \\textbf{%1.5f}" % (x))
        else:
            #sys.stdout.write(" & $%1.4f \\pm %1.4f$" % (x, y))
            sys.stdout.write(" & %1.5f" % (x))

    sys.stdout.write(" \\\\\n")

if __name__ == "__main__":
    mthds = ["CA-ML",
             "ASTER-wh",
             "TQMC-wh_n2",
             "TQMC-n2",
             "ASTRID-ws"]

    supps = ["abayes"]
    ngens = [50, 200, 1000]
    ntaxs = [10, 50, 100, 200, 500, 1000]

    sys.stdout.write("\\begin{table}[!h]\n")
    sys.stdout.write("\\caption[")
    sys.stdout.write("Speies tree error (RF) for ASTRAL-II simulated data]{\\textbf{Species tree (RF) error rate for ASTRAL-II simulated data.} ")
    sys.stdout.write("Mean error rate is given across replicates for each method. ")
    sys.stdout.write("All weighted methods were given gene trees with abayes support. ")
    sys.stdout.write("CA-ML trees from original study were not available for 1000-taxon, 1000-gene data sets.}\n")
    sys.stdout.write("\\centering\n")
    sys.stdout.write("\\scriptsize\n")
    sys.stdout.write("\\begin{tabular}{c c c c c c c c c c c}\n")
    sys.stdout.write("\\toprule \n")

    sys.stdout.write(" \\# of & ILS & speciation & \\# of & CA-ML & ASTER & TREE-QMC & TREE-QMC & ASTRID \\\\\n")
    sys.stdout.write(" taxa & level & location & genes &  & wh & wh (n2) & (n2) & ws \\\\\n")

    sys.stdout.write("\\midrule\n")
    for ntax in ntaxs:
        if ntax == 200:
            df = pandas.read_csv("data-varyils-error-and-qscore.csv", keep_default_na=True)
            hghts = ["low", "medium", "high"]
            rates = ["deep", "shallow"]
        else:
            df = pandas.read_csv("data-varyntax-error-and-qscore.csv", keep_default_na=True)
            hghts = ["medium"]
            rates = ["shallow"]

        for hght in hghts:
            for rate in rates:
                for ngen in ngens:
                    xdf = df[(df["NTAX"] == ntax) &
                             (df["ILSL"] == hght) &
                             (df["SPEC"] == rate) &
                             (df["NGEN"] == ngen)]

                    report_stats(xdf, mthds, ngen, ntax, hght, rate)


    sys.stdout.write("\\bottomrule\n")
    sys.stdout.write("\\end{tabular}\n")
    sys.stdout.write("\\end{table}\n")

"""
\begin{table}[!h]
\caption[Speies tree error (RF) for ASTRAL-II simulated data]{\textbf{Species tree (RF) error rate for ASTRAL-II simulated data.} Mean error rate is given across replicates for each method. All weighted methods were given gene trees with abayes support. CA-ML trees from original study were not available for 1000-taxon, 1000-gene data sets.}
\centering
\scriptsize
\begin{tabular}{c c c c c c c c c c c}
\toprule 
 \# of & ILS & speciation & \# of & CA-ML & ASTER & TREE-QMC & TREE-QMC & ASTRID \\
 taxa & level & location & genes &  & wh & wh (n2) & (n2) & ws \\
\midrule
10 & medium & shallow & 50 & 0.03830 & 0.03060 & 0.03060 & 0.03570 & \textbf{0.02810} \\
10 & medium & shallow & 200 & 0.01790 & \textbf{0.00760} & 0.01020 & 0.01790 & \textbf{0.00760} \\
10 & medium & shallow & 1000 & 0.02080 & \textbf{0.01560} & \textbf{0.01560} & \textbf{0.01560} & \textbf{0.01560} \\
50 & medium & shallow & 50 & 0.07800 & 0.06070 & 0.05850 & 0.07090 & \textbf{0.05760} \\
50 & medium & shallow & 200 & 0.04520 & 0.03060 & \textbf{0.02840} & 0.04120 & 0.03550 \\
50 & medium & shallow & 1000 & 0.02660 & 0.01860 & \textbf{0.01770} & 0.02530 & 0.01910 \\
100 & medium & shallow & 50 & 0.09120 & 0.07160 & \textbf{0.06720} & 0.07230 & 0.07040 \\
100 & medium & shallow & 200 & 0.04740 & 0.03830 & \textbf{0.03610} & 0.04680 & 0.03930 \\
100 & medium & shallow & 1000 & 0.02490 & \textbf{0.01910} & 0.01980 & 0.03040 & 0.02400 \\
200 & low & deep & 50 & \textbf{0.03980} & 0.05830 & 0.05160 & 0.06710 & 0.06490 \\
200 & low & deep & 200 & \textbf{0.02230} & 0.03520 & 0.03200 & 0.05060 & 0.05130 \\
200 & low & deep & 1000 & \textbf{0.01780} & 0.03010 & 0.02480 & 0.04520 & 0.04850 \\
200 & low & shallow & 50 & 0.05360 & 0.04810 & \textbf{0.04240} & 0.05280 & 0.04430 \\
200 & low & shallow & 200 & 0.03110 & 0.02350 & 0.02240 & 0.02910 & \textbf{0.02230} \\
200 & low & shallow & 1000 & 0.01440 & 0.01260 & \textbf{0.01120} & 0.01880 & 0.01300 \\
200 & medium & deep & 50 & 0.10300 & 0.08920 & \textbf{0.08570} & 0.09750 & 0.09100 \\
200 & medium & deep & 200 & 0.05690 & 0.05230 & \textbf{0.04660} & 0.05700 & 0.05500 \\
200 & medium & deep & 1000 & \textbf{0.02840} & 0.03520 & 0.03090 & 0.03980 & 0.03920 \\
200 & medium & shallow & 50 & 0.09220 & 0.07060 & \textbf{0.06690} & 0.08090 & 0.07280 \\
200 & medium & shallow & 200 & 0.05520 & 0.04030 & \textbf{0.03820} & 0.04750 & 0.04170 \\
200 & medium & shallow & 1000 & 0.02780 & 0.02410 & \textbf{0.02270} & 0.03170 & 0.02600 \\
200 & high & deep & 50 & 0.28160 & 0.18980 & \textbf{0.18580} & 0.20670 & 0.20630 \\
200 & high & deep & 200 & 0.16100 & \textbf{0.09130} & 0.09190 & 0.10940 & 0.10560 \\
200 & high & deep & 1000 & 0.08010 & \textbf{0.04710} & 0.04860 & 0.06170 & 0.04930 \\
200 & high & shallow & 50 & 0.27950 & 0.18180 & \textbf{0.17430} & 0.18980 & 0.20220 \\
200 & high & shallow & 200 & 0.16250 & 0.08940 & \textbf{0.08300} & 0.09600 & 0.10130 \\
200 & high & shallow & 1000 & 0.07950 & 0.04000 & \textbf{0.03770} & 0.04930 & 0.04710 \\
500 & medium & shallow & 50 & 0.09240 & 0.07400 & \textbf{0.06650} & 0.07630 & 0.07400 \\
500 & medium & shallow & 200 & 0.04720 & 0.04000 & \textbf{0.03410} & 0.04260 & 0.04100 \\
500 & medium & shallow & 1000 & 0.02330 & 0.02430 & \textbf{0.02030} & 0.02810 & 0.02610 \\
1000 & medium & shallow & 50 & 0.09760 & 0.08670 & \textbf{0.07680} & 0.10110 & 0.08810 \\
1000 & medium & shallow & 200 & 0.05150 & 0.04850 & \textbf{0.04180} & 0.05880 & 0.05180 \\
1000 & medium & shallow & 1000 & nan & 0.03050 & \textbf{0.02460} & 0.03640 & 0.03320 \\
\bottomrule
\end{tabular}
\end{table}
"""