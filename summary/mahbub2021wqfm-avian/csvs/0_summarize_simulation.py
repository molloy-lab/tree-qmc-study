import pandas
import numpy
import sys



def report_stats():
    """
    0.5X 500 : ILS = 0.590622 +/- 0.001831
    1X 500 : ILS = 0.473024 +/- 0.002017
    2X 500 : ILS = 0.353673 +/- 0.001791

    0.5X 500 : GTEE = 0.597427 +/- 0.093480
    1X 500 : GTEE = 0.597224 +/- 0.065515
    2X 500 : GTEE = 0.618598 +/- 0.002033

    0.5X 500 : AD = 0.687187 +/- 0.001360
    1X 500 : AD = 0.637170 +/- 0.001506
    2X 500 : AD = 0.596780 +/- 0.002313
    """
    scals = ["0.5X", "1X", "2X"]
    ngen = 1000
    nbps = 500
    repls = [str("R%d" % x) for x in range(1, 21)]

    # ILS
    df = pandas.read_csv("true_species_tree_vs_true_gene_trees.csv")
    for scal in scals:
        rfs = []
        for repl in repls:
            xdf = df[(df["SCAL"] == scal) &
                     (df["NGEN"] == ngen) &
                     (df["NBPS"] == nbps) &
                     (df["REPL"] == repl)]

            rfs.append(numpy.mean(xdf.RF.values))
        x = numpy.mean(rfs)
        y = numpy.std(rfs)
        print("%s %s : ILS = %f +/- %f" % (scal, nbps, x, y))

    # GTEE
    df = pandas.read_csv("true_vs_estimated_gene_trees.csv")
    for scal in scals:
        rfs = []
        for repl in repls:
            xdf = df[(df["SCAL"] == scal) &
                     (df["NGEN"] == ngen) &
                     (df["NBPS"] == nbps) &
                     (df["REPL"] == repl)]

            rfs.append(numpy.mean(xdf.RF.values))
        x = numpy.mean(rfs)
        y = numpy.std(rfs)
        print("%s %s : GTEE = %f +/- %f" % (scal, nbps, x, y))

    # AD
    df = pandas.read_csv("true_species_tree_vs_estimated_gene_trees.csv")
    for scal in scals:
        rfs = []
        for repl in repls:
            xdf = df[(df["SCAL"] == scal) &
                     (df["NGEN"] == ngen) &
                     (df["NBPS"] == nbps) &
                     (df["REPL"] == repl)]

            rfs.append(numpy.mean(xdf.RF.values))
        x = numpy.mean(rfs)
        y = numpy.std(rfs)
        print("%s %s : AD = %f +/- %f" % (scal, nbps, x, y))


if __name__ == "__main__":
    report_stats()
