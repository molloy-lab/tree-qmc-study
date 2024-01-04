import pandas
import numpy
import sys

sys.exit("DONE RUNNING")

def report_stats(df):
    """
    """
    nbpss = [200, 400, 800, 1600]
    repls = range(1, 51)

    for nbps in nbpss:
        keep = []

        for repl in repls:
            xdf = df[(df["NBPS"] == nbps) &
                     (df["REPL"] == repl)]

            fnrs = xdf.FN.values / xdf.I1.values

            keep.append(numpy.mean(fnrs))
    
        x = numpy.mean(keep)
        y = numpy.std(keep)
        print("%d : MEAN FNR = %f +/- %f" % (nbps, x, y))


if __name__ == "__main__":
    # ILS
    print("Evaluating ILS")
    df = pandas.read_csv("all_true_species_tree_vs_true_gene_trees.csv")
    report_stats(df)

    # GTEE
    print("Evaluating GTEE")
    df = pandas.read_csv("all_true_vs_estimated_gene_trees.csv")
    report_stats(df)

    # AD
    print("Evaluating AD")
    df = pandas.read_csv("all_true_species_tree_vs_estimated_gene_trees.csv")
    report_stats(df)

"""
Evaluating ILS
200 : MEAN FNR = 0.457506 +/- 0.060119
400 : MEAN FNR = 0.457506 +/- 0.060119
800 : MEAN FNR = 0.457506 +/- 0.060119
1600 : MEAN FNR = 0.457506 +/- 0.060119

Evaluating GTEE
200 : MEAN FNR = 0.554600 +/- 0.079195
400 : MEAN FNR = 0.422867 +/- 0.082320
800 : MEAN FNR = 0.311485 +/- 0.078926
1600 : MEAN FNR = 0.225808 +/- 0.073496

Evaluating AD
200 : MEAN FNR = 0.663386 +/- 0.065200
400 : MEAN FNR = 0.589428 +/- 0.065303
800 : MEAN FNR = 0.538503 +/- 0.062787
1600 : MEAN FNR = 0.507028 +/- 0.061079
"""