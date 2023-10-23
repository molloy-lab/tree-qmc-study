import pandas
import numpy
import sys

def report_stats():
    """
    noscale 250b : ILS = 0.329823 +/- 0.004795
    noscale 500b : ILS = 0.329823 +/- 0.004795
    noscale 1000b : ILS = 0.329823 +/- 0.004795
    noscale 1500b : ILS = 0.329823 +/- 0.004795
    
    noscale 250b : GTEE = 0.427331 +/- 0.009308
    noscale 500b : GTEE = 0.282221 +/- 0.057845
    noscale 1000b : GTEE = 0.160537 +/- 0.005587
    noscale 1500b : GTEE = 0.119206 +/- 0.005157

    noscale 250b : AD = 0.538199 +/- 0.006897
    noscale 500b : AD = 0.438250 +/- 0.005733
    noscale 1000b : AD = 0.379963 +/- 0.006049
    noscale 1500b : AD = 0.362118 +/- 0.005744
    """
    scal = "noscale"
    nbpss = ["250b", "500b", "1000b", "1500b"]
    repls = [str("R%d" % x) for x in range(1, 21)]

    # ILS
    df = pandas.read_csv("true_species_tree_vs_true_gene_trees.csv")
    for nbps in nbpss:
        rfs = []
        for repl in repls:
            xdf = df[(df["SCAL"] == scal) &
                     (df["NGEN"] == "200g") &
                     (df["NBPS"] == nbps) &
                     (df["REPL"] == repl)]
            rfs.append(numpy.mean(xdf.RF.values))
        x = numpy.mean(rfs)
        y = numpy.std(rfs)
        print("%s %s : ILS = %f +/- %f" % (scal, nbps, x, y))

    # GTEE
    df = pandas.read_csv("true_vs_estimated_gene_trees.csv")
    for nbps in nbpss:
        rfs = []
        for repl in repls:
            xdf = df[(df["SCAL"] == scal) &
                     (df["NGEN"] == "200g") &
                     (df["NBPS"] == nbps) &
                     (df["REPL"] == repl)]
            rfs.append(numpy.mean(xdf.RF.values))
        x = numpy.mean(rfs)
        y = numpy.std(rfs)
        print("%s %s : GTEE = %f +/- %f" % (scal, nbps, x, y))

    # AD
    df = pandas.read_csv("true_species_tree_vs_estimated_gene_trees.csv")
    for nbps in nbpss:
        rfs = []
        for repl in repls:
            xdf = df[(df["SCAL"] == scal) &
                     (df["NGEN"] == "200g") &
                     (df["NBPS"] == nbps) &
                     (df["REPL"] == repl)]
            rfs.append(numpy.mean(xdf.RF.values))
        x = numpy.mean(rfs)
        y = numpy.std(rfs)
        print("%s %s : AD = %f +/- %f" % (scal, nbps, x, y))


if __name__ == "__main__":
    report_stats()
