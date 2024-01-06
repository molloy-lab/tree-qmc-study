import pandas
import numpy
import sys

#sys.exit("DONE RUNNING")

def report_stats(df, mthds, supp, ngen, nbps):
    sys.stdout.write("%d & %d" % (ngen, nbps))

    xdf = df[(df["NGEN"] == ngen) &
             (df["NBPS"] == nbps) &
             (df["SUPP"] == supp)]

    ydf = df[(df["NGEN"] == ngen) &
             (df["NBPS"] == nbps) &
             (df["SUPP"] == "none_refinepoly")]  # Like original

    keep = []
    minval = 1
    for ind, mthd in enumerate(mthds):
        if mthd == "TQMC-wn_n2":
            rfs = ydf[(ydf["MTHD"] == mthd)].SERF.values
        else:
            rfs = xdf[(xdf["MTHD"] == mthd)].SERF.values

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
    sys.stdout.write("\\caption[RF error for S100 simulated data]{\\textbf{RF error rate for S100 simulated data.} ")
    sys.stdout.write("Mean error rate is given across 50 replicates for each method. ")
    sys.stdout.write("TREE-QMC-n2 does not use support; it is run with the refined abayes trees because the old TREE-QMC method could not handle polytomies}.")
    #sys.stdout.write("\\label{tab:}\n")
    sys.stdout.write("\\centering\n")
    sys.stdout.write("\\footnotesize\n")  # \small
    sys.stdout.write("\\begin{tabular}{r r c c c c c c c}\n")
    sys.stdout.write("\\toprule \n")

    sys.stdout.write(" \\# of & sequence & ASTRID & ASTER & TREE-QMC & TREE-QMC & TREE-QMC & TREE-QMC & TREE-QMC \\\\\n")
    sys.stdout.write(" genes & length & ws & wh & wh-n2 & wh-n1 & wh-n0 & ws-n2 & n2 \\\\\n")
    
    sys.stdout.write("\\midrule\n")

    mthds = ["ASTRID-ws",
             "ASTER-wh",
             "TQMC-wh_n2",
             "TQMC-wh_n1",
             "TQMC-wh_n0",
             "TQMC-ws_n2",
             "TQMC-wn_n2"]

    df = pandas.read_csv("data-all-error.csv", na_values='NA', keep_default_na=False)

    supps = ["bs", "abayes"]
    ngens = [50, 200, 500, 1000]
    nbpss = [200, 400, 800, 1600]

    for supp in supps:
        sys.stdout.write("\\multicolumn{5}{l}{\\textit{%s support for weighted methods}} \\\\[0.25em]\n" % supp)
        for ngen in ngens:
            for nbps in nbpss:
                report_stats(df, mthds, supp, ngen, nbps)

    sys.stdout.write("\\bottomrule\n")
    sys.stdout.write("\\end{tabular}\n")
    sys.stdout.write("\\end{table}\n")
