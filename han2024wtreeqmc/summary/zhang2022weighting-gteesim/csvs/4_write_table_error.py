import pandas
import numpy
import sys

#sys.exit("DONE RUNNING")

def report_stats(df, mthds, supp, ngen, nbps):
    sys.stdout.write("%d & %d" % (ngen, nbps))

    xdf = df[(df["NGEN"] == ngen) & (df["NBPS"] == nbps)]

    keep = []
    minval = 1
    for ind, mthd in enumerate(mthds):
        if mthd == "TQMC-n2":
            # IMPORTANT: TQMC-n2 doesn't use support
            # and bs support means unrefined trees
            # that are internally refined by TREE-QMC
            rfs = xdf[(xdf["MTHD"] == mthd) & (xdf["SUPP"] == "bs")].SERF.values
        else:
            rfs = xdf[(xdf["MTHD"] == mthd) & (xdf["SUPP"] == supp)].SERF.values

        if len(rfs) != 50:
            sys.exit("\nERROR")

        x = round(round(numpy.mean(rfs), 5), 4)
        #y = round(round(numpy.std(rfs), 5), 4)
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
    sys.stdout.write("\\begin{table}[!h]\n")
    sys.stdout.write("\\caption[Species tree error (RF) for S100 simulated data]{\\textbf{Species tree (RF) error rate for S100 simulated data.} ")
    sys.stdout.write("Mean error rate is given across 50 replicates for each method. ")
    sys.stdout.write("Note: TREE-QMC-n2 refines polytomies in the input gene trees randomly}.")
    #sys.stdout.write("\\label{tab:}\n")
    sys.stdout.write("\\centering\n")
    sys.stdout.write("\\footnotesize\n")  # \small
    sys.stdout.write("\\begin{tabular}{r r c c c c c c c}\n")
    sys.stdout.write("\\toprule \n")

    sys.stdout.write(" \\# of & sequence & ASTRID & ASTER & TREE-QMC & TREE-QMC & TREE-QMC & TREE-QMC & TREE-QMC \\\\\n")
    sys.stdout.write(" genes & length & ws & wh & wh-n2 & wh-n1 & wh-n0 & ws-n2 & n2 \\\\\n")

    mthds = ["ASTRID-ws",
             "ASTER-wh",
             "TQMC-wh_n2",
             "TQMC-wh_n1",
             "TQMC-wh_n0",
             "TQMC-ws_n2",
             "TQMC-n2"]

    df = pandas.read_csv("data-all-error-and-qscore.csv", na_values='NA', keep_default_na=False)

    supps = ["bs", "abayes"]
    ngens = [50, 200, 500, 1000]
    nbpss = [200, 400, 800, 1600]

    for supp in supps:
        sys.stdout.write("\\midrule\n")
        sys.stdout.write("\\multicolumn{9}{c}{\\textbf{%s support for weighted methods}} \\\\\n" % supp)
        sys.stdout.write("\\midrule\n")
        for ngen in ngens:
            for nbps in nbpss:
                report_stats(df, mthds, supp, ngen, nbps)

    sys.stdout.write("\\bottomrule\n")
    sys.stdout.write("\\end{tabular}\n")
    sys.stdout.write("\\end{table}\n")


"""
\begin{table}[!h]
\caption[Species tree error (RF) for S100 simulated data]{\textbf{Species tree (RF) error rate for S100 simulated data.} Mean error rate is given across 50 replicates for each method. Note: TREE-QMC-n2 refines polytomies in the input gene trees randomly}.\centering
\footnotesize
\begin{tabular}{r r c c c c c c c}
\toprule 
 \# of & sequence & ASTRID & ASTER & TREE-QMC & TREE-QMC & TREE-QMC & TREE-QMC & TREE-QMC \\
 genes & length & ws & wh & wh-n2 & wh-n1 & wh-n0 & ws-n2 & n2 \\
\midrule
\multicolumn{9}{c}{\textbf{bs support for weighted methods}} \\
\midrule
50 & 200 & 0.16310 & 0.15530 & \textbf{0.15390} & 0.15630 & 0.16740 & 0.15610 & 0.17590 \\
50 & 400 & 0.13120 & 0.12550 & \textbf{0.12220} & 0.12350 & 0.13200 & 0.12350 & 0.13100 \\
50 & 800 & 0.11350 & 0.10390 & \textbf{0.10310} & 0.10430 & 0.11040 & 0.10530 & 0.11450 \\
50 & 1600 & 0.10350 & 0.09590 & \textbf{0.09290} & 0.09530 & 0.09920 & 0.09690 & 0.09710 \\
200 & 200 & 0.09710 & 0.09630 & \textbf{0.09490} & 0.09880 & 0.10020 & 0.09820 & 0.10960 \\
200 & 400 & 0.07960 & 0.07670 & \textbf{0.07290} & 0.07690 & 0.08330 & 0.07470 & 0.08550 \\
200 & 800 & 0.06900 & 0.06470 & \textbf{0.06290} & 0.06430 & 0.06800 & 0.06410 & 0.06980 \\
200 & 1600 & 0.06140 & 0.05820 & \textbf{0.05550} & 0.05880 & 0.06160 & 0.05610 & 0.06410 \\
500 & 200 & 0.07740 & \textbf{0.07530} & 0.07550 & 0.07920 & 0.08020 & 0.07940 & 0.09260 \\
500 & 400 & 0.06160 & 0.05960 & \textbf{0.05550} & 0.05960 & 0.06450 & 0.06000 & 0.07220 \\
500 & 800 & 0.05080 & 0.04900 & \textbf{0.04800} & 0.04960 & 0.05350 & 0.04920 & 0.05470 \\
500 & 1600 & 0.04410 & 0.04290 & \textbf{0.04060} & 0.04260 & 0.04430 & 0.04240 & 0.04880 \\
1000 & 200 & 0.06690 & \textbf{0.06530} & 0.06730 & 0.07160 & 0.07020 & 0.07040 & 0.08000 \\
1000 & 400 & 0.05160 & 0.05140 & \textbf{0.04800} & 0.05020 & 0.05350 & 0.05060 & 0.06330 \\
1000 & 800 & 0.04200 & 0.04180 & \textbf{0.03900} & 0.04260 & 0.04490 & 0.04120 & 0.05040 \\
1000 & 1600 & 0.03590 & 0.03610 & \textbf{0.03510} & 0.03630 & 0.03820 & 0.03530 & 0.04100 \\
\midrule
\multicolumn{9}{c}{\textbf{abayes support for weighted methods}} \\
\midrule
50 & 200 & 0.15040 & 0.14570 & \textbf{0.14140} & 0.14570 & 0.15860 & 0.14550 & 0.17590 \\
50 & 400 & 0.12390 & 0.12180 & \textbf{0.11120} & 0.11690 & 0.12330 & 0.11530 & 0.13100 \\
50 & 800 & 0.11140 & 0.10840 & \textbf{0.09670} & 0.10370 & 0.11260 & 0.10060 & 0.11450 \\
50 & 1600 & 0.10450 & 0.09550 & \textbf{0.08630} & 0.09310 & 0.09940 & 0.08820 & 0.09710 \\
200 & 200 & 0.09470 & 0.09220 & \textbf{0.08310} & 0.08710 & 0.09800 & 0.09430 & 0.10960 \\
200 & 400 & 0.07430 & 0.07220 & \textbf{0.06710} & 0.07120 & 0.07780 & 0.07180 & 0.08550 \\
200 & 800 & 0.06780 & 0.06410 & \textbf{0.06180} & 0.06410 & 0.06760 & 0.06240 & 0.06980 \\
200 & 1600 & 0.05880 & 0.05650 & \textbf{0.05160} & 0.05490 & 0.05800 & 0.05370 & 0.06410 \\
500 & 200 & 0.07530 & 0.07160 & \textbf{0.06920} & 0.07290 & 0.07780 & 0.07840 & 0.09260 \\
500 & 400 & 0.06080 & \textbf{0.05760} & 0.05860 & 0.05840 & 0.06080 & 0.06160 & 0.07220 \\
500 & 800 & 0.04880 & 0.04690 & \textbf{0.04370} & 0.04590 & 0.04960 & 0.04820 & 0.05470 \\
500 & 1600 & 0.04080 & 0.03920 & \textbf{0.03530} & 0.03860 & 0.04220 & 0.03820 & 0.04880 \\
1000 & 200 & 0.06220 & 0.06310 & \textbf{0.05350} & 0.06100 & 0.06840 & 0.05920 & 0.08000 \\
1000 & 400 & 0.05120 & 0.05160 & \textbf{0.04650} & 0.05160 & 0.05610 & 0.05240 & 0.06330 \\
1000 & 800 & 0.04220 & 0.04180 & \textbf{0.03630} & 0.04060 & 0.04650 & 0.04140 & 0.05040 \\
1000 & 1600 & 0.03740 & 0.03530 & \textbf{0.02880} & 0.03490 & 0.04000 & 0.03410 & 0.04100 \\
\bottomrule
\end{tabular}
\end{table}
"""