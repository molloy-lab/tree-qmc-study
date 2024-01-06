import pandas
import numpy
import sys

sys.exit("DONE RUNNING")

def report_stats(df, mthds, ntax, ngen, nbps, blsc, psiz, miss, vary, emet):
    sys.stdout.write("%s" % vary)

    xdf = df[(df["NTAX"] == ntax) &
             (df["NGEN"] == ngen) &
             (df["NBPS"] == nbps) &
             (df["BLSC"] == blsc) &
             (df["PSIZ"] == psiz) &
             (df["MISS"] == miss)]

    for mthd in mthds:
        if emet == "FNR":
            er = xdf[(xdf["MTHD"] == mthd)].SEFNR.values
        elif emet == "FPR":
            er = xdf[(xdf["MTHD"] == mthd)].SEFPR.values
        else:
            sys.exit("unknown error metric\n")

        if len(er) != 50:
            sys.exit("ERROR")

        x = round(round(numpy.mean(er), 5), 4)
        y = round(round(numpy.std(er), 5), 4)
        #sys.stdout.write(" & $%1.4f \\pm %1.4f$" % (x, y))
        sys.stdout.write(" & %1.4f" % (x))

    sys.stdout.write(" \\\\\n")

if __name__ == "__main__":
    
    emet = "FNR"
    emet = "FPR"

    sys.stdout.write("\\begin{table}[!h]\n")
    sys.stdout.write("\\caption[")
    if emet == "FNR":
        sys.stdout.write("FNR for Asteroid simulated data]{\\textbf{False negative error rate for Asteroid simulated data.} ")
    else:
        sys.stdout.write("FPR for Asteroid simulated data]{\\textbf{False positive error rate for Asteroid simulated data.} ")
    sys.stdout.write("Mean error rate is given across 50 replicates for each method.}")
    #sys.stdout.write("\\label{tab:}\n")
    sys.stdout.write("\\centering\n")
    sys.stdout.write("\\small\n")
    sys.stdout.write("\\begin{tabular}{r c c c c c c c c}\n")
    sys.stdout.write("\\toprule \n")

    sys.stdout.write(" & ASTRID & ASTER & Asteroid & TREE-QMC-n2 & n1 & n0 & n2-shared \\\\\n")
    sys.stdout.write("\\midrule\n")

    experiments = ["varypsiz", "varyntax", "varyngen",
                   "varynbps", "varyblsc", "varymiss"]

    mthds = ["wastrid_vanilla",
             "aster_v1.16.3.4",
             "asteroid",
             "wtreeqmc_wf_n2",
             "wtreeqmc_wf_n1",
             "wtreeqmc_wf_n0",
             "wtreeqmc_wf_n2_shared"]

    for do in experiments:
        ntax = 50
        ngen = 1000
        nbps = 100
        blsc = 1.0
        psiz = 50000000
        miss = 0.6  # ms always same as mf

        df = pandas.read_csv("data-" + do + "-error-and-timings.csv",
                             na_values='NA', keep_default_na=False)

        if do == "varypsiz":
            sys.stdout.write("\\multicolumn{5}{l}{\\textit{Varying population size}} \\\\[0.25em]\n")
            names = ["        10", 
                     "  50000000",
                     " 100000000",
                     " 500000000",
                     "1000000000"]
            psizs = [10, 50000000, 100000000, 500000000, 1000000000]
            for name, psiz in zip(names, psizs):
                report_stats(df, mthds, ntax, ngen, nbps, blsc, psiz, miss, name, emet)
        elif do == "varyntax":
            sys.stdout.write("\\multicolumn{5}{l}{\\textit{Varying number of taxa}} \\\\[0.25em]\n")
            names = ["  25", "  75", "  50", " 100", " 125", " 150"]
            ntaxs = [25, 75, 50, 100, 125, 150]  # dup 50
            for name, ntax in zip(names, ntaxs):
                report_stats(df, mthds, ntax, ngen, nbps, blsc, psiz, miss, name, emet)
        elif do == "varyngen":
            sys.stdout.write("\\multicolumn{5}{l}{\\textit{Varying number of genes}} \\\\[0.25em]\n")
            names = [" 250", " 500", "1000", "2000"]
            ngens = [250, 500, 1000, 2000]  # dup 1000
            for name, ngen in zip(names, ngens):
                report_stats(df, mthds, ntax, ngen, nbps, blsc, psiz, miss, name, emet)
        elif do == "varynbps":
            sys.stdout.write("\\multicolumn{5}{l}{\\textit{Varying sequence length}} \\\\[0.25em]\n")
            names = ["  50", " 100", " 200", " 500"]
            nbpss = [50, 100, 200, 500]  # dup 100
            for name, nbps in zip(names, nbpss):
                report_stats(df, mthds, ntax, ngen, nbps, blsc, psiz, miss, name, emet)
        elif do == "varyblsc":
            sys.stdout.write("\\multicolumn{5}{l}{\\textit{Varying branch length scaler}} \\\\[0.25em]\n")
            names = ["  0.05", "  0.10", "  1.00", " 10.00", "100.00", "200.00"]
            blscs = [0.05, 0.1, 1, 10, 100, 200]  # dup 1.0
            for name, blsc in zip(names, blscs):
                report_stats(df, mthds, ntax, ngen, nbps, blsc, psiz, miss, name, emet)
        elif do == "varymiss":
            sys.stdout.write("\\multicolumn{5}{l}{\\textit{Varying missingness parameter}} \\\\[0.25em]\n")
            names = ["0.50", "0.55", "0.60", "0.65", "0.70", "0.75"]
            misss = [0.5, 0.55, 0.6, 0.65, 0.7, 0.75]  # dup 0.6
            for name, miss in zip(names, misss):
                report_stats(df, mthds, ntax, ngen, nbps, blsc, psiz, miss, name, emet)
        else:
            sys.exit()

    sys.stdout.write("\\bottomrule\n")
    sys.stdout.write("\\end{tabular}\n")
    sys.stdout.write("\\end{table}\n")

