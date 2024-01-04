import pandas
import numpy
import sys

sys.exit("DONE RUNNING")

def report_stats(df):
    """
    """
    repls = range(1, 51)

    keep = []

    repls = list(set(df.REPL.values))
    for repl in repls:
        xdf = df[df["REPL"] == repl]
        fnrs = xdf.FN.values / xdf.I1.values
        keep.append(numpy.mean(fnrs))
    x = numpy.mean(keep)
    y = numpy.std(keep)
    sys.stdout.write(" & %1.4f +/- %1.4f" % (x, y))


if __name__ == "__main__":
    nbpss = [200, 400, 800, 1600]

    df_ils = pandas.read_csv("all_true_species_tree_vs_true_gene_trees.csv")
    df_gtee = pandas.read_csv("all_true_vs_estimated_gene_trees.csv")
    df_ad = pandas.read_csv("all_true_species_tree_vs_estimated_gene_trees.csv")

    sys.stdout.write("NBPS & ILS & GTEE & AD\\\\\n")
    for nbps in nbpss:
        sys.stdout.write("%d" % nbps)

        xdf_ils = df_ils[df_ils["NBPS"] == nbps]
        report_stats(xdf_ils)
        
        xdf_gtee = df_gtee[df_gtee["NBPS"] == nbps]
        report_stats(xdf_gtee)

        xdf_ad = df_ad[df_ad["NBPS"] == nbps]
        report_stats(xdf_ad)

        sys.stdout.write("\n")


"""
NBPS & ILS & GTEE & AD\\
200 & 0.4575 +/- 0.0601 & 0.5546 +/- 0.0792 & 0.6634 +/- 0.0652
400 & 0.4575 +/- 0.0601 & 0.4229 +/- 0.0823 & 0.5894 +/- 0.0653
800 & 0.4575 +/- 0.0601 & 0.3115 +/- 0.0789 & 0.5385 +/- 0.0628
1600 & 0.4575 +/- 0.0601 & 0.2258 +/- 0.0735 & 0.5070 +/- 0.0611
"""
