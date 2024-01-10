import pandas
import numpy
import sys

#sys.exit("DONE RUNNING")

def report_stats(df, mthds, ntax, ngen, nbps, blsc, psiz, miss, vary, emet, doit):
    sys.stdout.write("& %s" % vary)

    xdf = df[(df["NTAX"] == ntax) &
             (df["NGEN"] == ngen) &
             (df["NBPS"] == nbps) &
             (df["BLSC"] == blsc) &
             (df["PSIZ"] == psiz) &
             (df["MISS"] == miss)]

    keep = []
    minval = 1
    for ind, mthd in enumerate(mthds):
        if emet == "FNR":
            ers = xdf[(xdf["MTHD"] == mthd)].SEFNR.values
        elif emet == "FPR":
            ers = xdf[(xdf["MTHD"] == mthd)].SEFPR.values
        else:
            sys.exit("unknown error metric\n")

        if len(ers) != 50:
            sys.exit("ERROR")

        x = round(round(numpy.mean(ers), 5), 4)
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

    if ((doit != "varypsiz") and \
        (ntax == 50) and \
        (ngen == 1000) and \
        (nbps == 100) and \
        (blsc == 1.0) and \
        (psiz == 50000000) and \
        (miss == 0.6)):
        sys.stdout.write(" (dup)")

    sys.stdout.write(" \\\\\n")


if __name__ == "__main__":
    
    emet = "FNR"
    emet = "FPR"

    sys.stdout.write("\\begin{table}[!h]\n")
    sys.stdout.write("\\caption[")
    if emet == "FNR":
        sys.stdout.write("Species tree error (FN) for Asteroid data]{\\textbf{Species tree (FN) error rate for Asteroid data.} ")
    else:
        sys.stdout.write("Species tree error (FP) for Asteroid data]{\\textbf{Species tree (FP) error rate for Asteroid data.} ")
    sys.stdout.write("Mean error rate is given across 50 replicates for each method.}\n")
    #sys.stdout.write("\\label{tab:}\n")
    sys.stdout.write("\\centering\n")
    sys.stdout.write("\\footnotesize\n")
    #sys.stdout.write("\\begin{tabular}{r c c c c c c c l}\n")
    sys.stdout.write("\\begin{tabular}{r c c c c c c c l}\n")
    sys.stdout.write("\\toprule \n")

    sys.stdout.write(" & & ASTRID & ASTER & Asteroid & TREE-QMC-n2 & n1 & n0 & n2-shared \\\\\n")
    sys.stdout.write("\\midrule\n")

    experiments = ["varypsiz", "varyntax", "varyngen",
                   "varynbps", "varyblsc", "varymiss"]

    mthds = ["ASTRID",
             "ASTER",
             "Asteroid",
             "TREE-QMC-wf_n2",
             "TREE-QMC-wf_n1",
             "TREE-QMC-wf_n0",
             "TREE-QMC-wf_n1_shared"]

    for doit in experiments:
        sys.stdout.write("\\midrule\n")

        ntax = 50
        ngen = 1000
        nbps = 100
        blsc = 1.0
        psiz = 50000000
        miss = 0.6  # ms always same as mf

        df = pandas.read_csv("data-" + doit + "-error-and-timings.csv",
                             na_values='NA', keep_default_na=False)

        if doit == "varypsiz":
            #sys.stdout.write("\\multicolumn{5}{l}{\\textit{Varying population size}} \\\\[0.25em]\n")
            sys.stdout.write("\\multirow{ 5}{2cm}{Varying population size}\n")
            names = ["        10", 
                     "  50000000",
                     " 100000000",
                     " 500000000",
                     "1000000000"]
            psizs = [10, 50000000, 100000000, 500000000, 1000000000]
            for name, psiz in zip(names, psizs):
                report_stats(df, mthds, ntax, ngen, nbps, blsc, psiz, miss, name, emet, doit)
        elif doit == "varyntax":
            #sys.stdout.write("\\multicolumn{5}{l}{\\textit{Varying number of taxa}} \\\\[0.25em]\n")
            sys.stdout.write("\\multirow{ 6}{2cm}{Varying number of taxa}\n")
            names = ["  25", "  75", "  50", " 100", " 125", " 150"]
            ntaxs = [25, 75, 50, 100, 125, 150]  # dup 50
            for name, ntax in zip(names, ntaxs):
                report_stats(df, mthds, ntax, ngen, nbps, blsc, psiz, miss, name, emet, doit)
        elif doit == "varyngen":
            #sys.stdout.write("\\multicolumn{5}{l}{\\textit{Varying number of genes}} \\\\[0.25em]\n")
            sys.stdout.write("\\multirow{ 4}{2cm}{Varying number of genes}\n")
            names = [" 250", " 500", "1000", "2000"]
            ngens = [250, 500, 1000, 2000]  # dup 1000
            for name, ngen in zip(names, ngens):
                report_stats(df, mthds, ntax, ngen, nbps, blsc, psiz, miss, name, emet, doit)
        elif doit == "varynbps":
            #sys.stdout.write("\\multicolumn{5}{l}{\\textit{Varying sequence length}} \\\\[0.25em]\n")
            sys.stdout.write("\\multirow{ 4}{2cm}{Varying sequence length}\n")
            names = ["  50", " 100", " 200", " 500"]
            nbpss = [50, 100, 200, 500]  # dup 100
            for name, nbps in zip(names, nbpss):
                report_stats(df, mthds, ntax, ngen, nbps, blsc, psiz, miss, name, emet, doit)
        elif doit == "varyblsc":
            #sys.stdout.write("\\multicolumn{5}{l}{\\textit{Varying branch length scaler}} \\\\[0.25em]\n")
            sys.stdout.write("\\multirow{ 6}{2cm}{Varying branch length scaler}\n")
            names = ["  0.05", "  0.10", "  1.00", " 10.00", "100.00", "200.00"]
            blscs = [0.05, 0.1, 1, 10, 100, 200]  # dup 1.0
            for name, blsc in zip(names, blscs):
                report_stats(df, mthds, ntax, ngen, nbps, blsc, psiz, miss, name, emet, doit)
        elif doit == "varymiss":
            #sys.stdout.write("\\multicolumn{5}{l}{\\textit{Varying missingness parameter}} \\\\[0.25em]\n")
            sys.stdout.write("\\multirow{ 6}{2cm}{Varying missingness parameter}\n")
            names = ["0.50", "0.55", "0.60", "0.65", "0.70", "0.75"]
            misss = [0.5, 0.55, 0.6, 0.65, 0.7, 0.75]  # dup 0.6
            for name, miss in zip(names, misss):
                report_stats(df, mthds, ntax, ngen, nbps, blsc, psiz, miss, name, emet, doit)
        else:
            sys.exit()

    sys.stdout.write("\\bottomrule\n")
    sys.stdout.write("\\end{tabular}\n")
    sys.stdout.write("\\end{table}\n")

