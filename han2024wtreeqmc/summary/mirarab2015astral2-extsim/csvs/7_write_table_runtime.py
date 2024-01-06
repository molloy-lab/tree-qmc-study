import pandas
import numpy
import sys
import warnings
import csv

def report_stats(df, mthds, ngen, ntax, hght, rate, supp):
    sys.stdout.write("%d & %s & %s & %d" % (ntax, hght, rate, ngen))

    nrepl = len(df[(df["MTHD"] == "TQMC-wh_n2")].REPL.values)

    for mthd in mthds:
        xdf = df[(df["MTHD"] == mthd) & (df["SUPP"] == supp)]

        vals = xdf.SECS.values

        if len(vals) != nrepl:
            sys.exit("\nERROR - inconsistent number of replicates")
        if len(vals) > 50:
            sys.exit("\nERROR - too many replicates")

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            fnravg = numpy.nanmean(vals)

        x = round(round(fnravg, 2), 1)
        sys.stdout.write(" & %1.1f" % (x))

    sys.stdout.write(" \\\\\n")

if __name__ == "__main__":
    mthds = ["ASTRID-ws",
             "TQMC-wh_n2",
             "ASTER-wh (1 thread)",
             "ASTER-wh (16 threads)"]

    supp = "sh"

    sys.stdout.write("\\begin{table}[!h]\n")
    sys.stdout.write("\\caption[")
    sys.stdout.write("Runtime for ASTRAL-II simulated data]{\\textbf{Runtime for ASTRAL-II simulated data.} ")
    sys.stdout.write("Mean runtime (in seconds) is given across replicates for each method. ")
    sys.stdout.write("All weighted methods were given gene trees with \\textbf{%s} support.}" % supp)
    #sys.stdout.write("\\label{tab:}\n")
    sys.stdout.write("\\centering\n")
    sys.stdout.write("\\scriptsize\n")
    sys.stdout.write("\\begin{tabular}{c c c c r r r r}\n")
    sys.stdout.write("\\toprule \n")

    sys.stdout.write(" \\# of & ILS & speciation & \\# of & ASTRID & TREE-QMC & ASTER & ASTER \\\\\n")
    sys.stdout.write(" taxa & level & location & genes & ws & wh-n2 & wh & wh  \\\\\n")
    sys.stdout.write("  &  &  &  &  &  & (1 thread) & (16 threads)  \\\\\n")
    sys.stdout.write("\\midrule\n")

    
    ngens = [1000]
    ntaxs = [10, 50, 100, 200, 500, 1000]
    for ntax in ntaxs:
        if ntax == 200:
            df = pandas.read_csv("data-varyils-error-and-runtime.csv",
                                 keep_default_na=True, quoting=csv.QUOTE_NONNUMERIC)
            hghts = ["low", "medium", "high"]
            rates = ["deep", "shallow"]
        else:
            df = pandas.read_csv("data-varyntax-error-and-runtime.csv",
                                 keep_default_na=True, quoting=csv.QUOTE_NONNUMERIC)
            hghts = ["medium"]
            rates = ["shallow"]

        for hght in hghts:
            for rate in rates:
                for ngen in ngens:
                    xdf = df[(df["NTAX"] == ntax) &
                             (df["ILSL"] == hght) &
                             (df["SPEC"] == rate) &
                             (df["NGEN"] == ngen)]
                    report_stats(xdf, mthds, ngen, ntax, hght, rate, supp)

    sys.stdout.write("\\bottomrule\n")
    sys.stdout.write("\\end{tabular}\n")
    sys.stdout.write("\\end{table}\n")

