import pandas
import numpy
import sys
import warnings
import csv

def report_stats(df, mthds, ngen, ntax, hght, rate):
    sys.stdout.write("%d & %s & %s & %d" % (ntax, hght, rate, ngen))

    nrepl = len(df[(df["MTHD"] == "TQMC-wh_n2")].REPL.values)

    for mthd in mthds:
        xdf = df[(df["MTHD"] == mthd) & (df["SUPP"] == "abayes")]

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
             "TQMC-n2",
             "TQMC-wh_n2",
             "ASTER-wh (16 threads)",
             "ASTER-wh (1 thread)"]

    sys.stdout.write("\\begin{table}[!h]\n")
    sys.stdout.write("\\caption[")
    sys.stdout.write("Runtime for ASTRAL-II simulated data]{\\textbf{Runtime for ASTRAL-II simulated data.} ")
    sys.stdout.write("Mean runtime (in seconds) is given across replicates for each method. ")
    sys.stdout.write("All weighted methods were given gene trees with abayes support. " )
    sys.stdout.write("CA-ML trees from original study were not available for 1000-taxon, 1000-gene data sets.}\n")
    sys.stdout.write("\\centering\n")
    sys.stdout.write("\\scriptsize\n")
    sys.stdout.write("\\begin{tabular}{c c c c r r r r r}\n")
    sys.stdout.write("\\toprule \n")

    sys.stdout.write(" \\# of & ILS & speciation & \\# of & ASTRID & TREE-QMC & TREE-QMC & ASTER & ASTER \\\\\n")
    sys.stdout.write(" taxa & level & location & genes & ws & (n2) & wh (n2) & wh & wh  \\\\\n")
    sys.stdout.write("  &  &  & &  &  &  & (16 threads) & (1 thread)  \\\\\n")
    sys.stdout.write("\\midrule\n")
    
    ngens = [1000]
    ntaxs = [10, 50, 100, 200, 500, 1000]
    for ntax in ntaxs:
        if ntax == 200:
            df = pandas.read_csv("data-varyils-runtime.csv",
                                 keep_default_na=True, quoting=csv.QUOTE_NONNUMERIC)
            hghts = ["low", "medium", "high"]
            rates = ["deep", "shallow"]
        else:
            df = pandas.read_csv("data-varyntax-runtime.csv",
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
                    report_stats(xdf, mthds, ngen, ntax, hght, rate)

    sys.stdout.write("\\bottomrule\n")
    sys.stdout.write("\\end{tabular}\n")
    sys.stdout.write("\\end{table}\n")

"""
\begin{table}[!h]
\caption[Runtime for ASTRAL-II simulated data]{\textbf{Runtime for ASTRAL-II simulated data.} Mean runtime (in seconds) is given across replicates for each method. All weighted methods were given gene trees with abayes support. CA-ML trees from original study were not available for 1000-taxon, 1000-gene data sets.}
\centering
\scriptsize
\begin{tabular}{c c c c r r r r r}
\toprule 
 \# of & ILS & speciation & \# of & ASTRID & TREE-QMC & TREE-QMC & ASTER & ASTER \\
 taxa & level & location & genes & ws & (n2) & wh (n2) & wh & wh  \\
  &  &  & &  &  &  & (16 threads) & (1 thread)  \\
\midrule
10 & medium & shallow & 1000 & 0.0 & 0.2 & 0.3 & 0.1 & 0.7 \\
50 & medium & shallow & 1000 & 0.1 & 3.3 & 7.5 & 2.4 & 34.5 \\
100 & medium & shallow & 1000 & 0.1 & 13.8 & 40.5 & 10.7 & 131.6 \\
200 & low & deep & 1000 & 0.3 & 51.8 & 170.8 & 41.5 & 700.0 \\
200 & low & shallow & 1000 & 0.3 & 50.4 & 159.9 & 36.9 & 645.5 \\
200 & medium & deep & 1000 & 0.3 & 55.7 & 181.3 & 44.0 & 734.3 \\
200 & medium & shallow & 1000 & 0.3 & 51.7 & 167.1 & 41.1 & 661.1 \\
200 & high & deep & 1000 & 0.4 & 61.1 & 201.6 & 54.0 & 955.0 \\
200 & high & shallow & 1000 & 0.4 & 58.8 & 192.7 & 53.8 & 869.3 \\
500 & medium & shallow & 1000 & 1.9 & 322.8 & 1084.0 & 273.1 & 5145.7 \\
1000 & medium & shallow & 1000 & 7.3 & 1303.3 & 4469.8 & 1272.5 & 28758.8 \\
\bottomrule
\end{tabular}
\end{table}
"""