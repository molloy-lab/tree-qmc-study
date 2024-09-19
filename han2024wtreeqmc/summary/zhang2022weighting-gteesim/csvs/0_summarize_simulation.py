import pandas
import numpy
import sys

sys.exit("DONE RUNNING")

def report_stats(df):
    repls = range(1, 51)
    keep = []
    repls = list(set(df.REPL.values))
    for repl in repls:
        xdf = df[df["REPL"] == repl]
        fnrs = xdf.FN.values / xdf.I1.values
        keep.append(numpy.mean(fnrs))
    x = numpy.mean(keep)
    y = numpy.std(keep)
    sys.stdout.write(" & $%1.4f \\pm %1.4f$" % (x, y))

def report_stats_poly(df):
    repls = range(1, 51)
    keep = []
    repls = list(set(df.REPL.values))
    for repl in repls:
        xdf = df[df["REPL"] == repl]
        ncon = xdf.I1.values - xdf.I2.values # Number of contracted branches
        keep.append(numpy.mean(ncon))
    x = numpy.mean(keep)
    y = numpy.std(keep)
    sys.stdout.write(" & $%1.4f \\pm %1.4f$" % (x, y))

if __name__ == "__main__":
    nbpss = [200, 400, 800, 1600]

    df_ils = pandas.read_csv("all_true_species_tree_vs_true_gene_trees.csv.gz", compression="gzip")
    df_gtee = pandas.read_csv("all_true_vs_estimated_gene_trees.csv.gz", compression="gzip")
    df_ad = pandas.read_csv("all_true_species_tree_vs_estimated_gene_trees.csv.gz", compression="gzip")

    sys.stdout.write("\\begin{table}[!h]\n")
    sys.stdout.write("\\caption[Properties of S100 simulated data]")
    sys.stdout.write("{\\textbf{Properties of S100 simulated data from \\cite{zhang2018astral3}.} ")
    sys.stdout.write("Number of refinements is the number of refinements needed to make FastTrees binary (they can be non-binary due to identical sequences).}\n")
    sys.stdout.write("\\label{tab:mahbub}\n")
    sys.stdout.write("\\centering\n")
    sys.stdout.write("\\small\n")
    sys.stdout.write("\\begin{tabular}{r c c c c}\n")
    sys.stdout.write("\\toprule \n")

    sys.stdout.write("Sequence length & ILS & GTEE & AD & \\# of refinements\\\\\n")
    sys.stdout.write("\\midrule\n")

    for nbps in nbpss:
        sys.stdout.write("%d" % nbps)

        xdf_ils = df_ils[df_ils["NBPS"] == nbps]
        report_stats(xdf_ils)

        xdf_gtee = df_gtee[df_gtee["NBPS"] == nbps]
        report_stats(xdf_gtee)

        xdf_ad = df_ad[df_ad["NBPS"] == nbps]
        report_stats(xdf_ad)

        xdf_ad = df_ad[df_ad["NBPS"] == nbps]
        report_stats_poly(xdf_ad)

        sys.stdout.write(" \\\\\n")

    sys.stdout.write("\\bottomrule\n")
    sys.stdout.write("\\end{tabular}\n")
    sys.stdout.write("\\end{table}\n")


"""
\begin{table}[!h]
\caption[Properties of S100 simulated data]{\textbf{Properties of S100 simulated data from \cite{zhang2018astral3}.} Number of refinements is the number of refinements needed to make FastTrees binary (they can be non-binary due to identical sequences).}
\label{tab:mahbub}
\centering
\small
\begin{tabular}{r c c c c}
\toprule 
Sequence length & ILS & GTEE & AD & \# of refinements\\
\midrule
200 & $0.4575 \pm 0.0601$ & $0.5546 \pm 0.0792$ & $0.6634 \pm 0.0652$ & $2.6337 \pm 3.9622$ \\
400 & $0.4575 \pm 0.0601$ & $0.4229 \pm 0.0823$ & $0.5894 \pm 0.0653$ & $0.8408 \pm 1.5477$ \\
800 & $0.4575 \pm 0.0601$ & $0.3115 \pm 0.0789$ & $0.5385 \pm 0.0628$ & $0.2706 \pm 0.6069$ \\
1600 & $0.4575 \pm 0.0601$ & $0.2258 \pm 0.0735$ & $0.5070 \pm 0.0611$ & $0.0972 \pm 0.2746$ \\
\bottomrule
\end{tabular}
\end{table}
"""
