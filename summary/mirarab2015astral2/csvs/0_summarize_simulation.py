import pandas
import numpy
import sys

"""
From ASTRAL-II paper:
FastTree can output polytomies when sequence alignments cannot 
distinguish between competing tree resolutions. We removed any 
gene tree where more than 50% of the internal nodes were polytomies. 
This pruning left fewer than 500 genes for three replicates of the 
200 taxon/500 K/1e-06 and 50-taxon model conditions, two replicates 
of the 100-taxon model condition and one replicate of the 10-taxon 
model condition. Those nine replicates (out of 550) were removed 
from our analyses.

To summarize:
+ remove 1 replicate (i.e. 41) from 10-taxon
+ remove 3 replicates (i.e. 21 and 41) from 50-taxon (next is 27 with 590ish genes)
+ remove 2 replicates (i.e. 8 and 47) from 100-taxon
+ remove 3 replicates (i.e. 8, 15, 49) from 200 taxon/500 K/1e-06
"""

def find_bad_repls_200tax():
    """
    ntax = 200, strht = 1e+07, srate = 1e-07
    ntax = 200, strht = 1e+07, srate = 1e-06
    ntax = 200, strht = 2e+06, srate = 1e-07
    ntax = 200, strht = 2e+06, srate = 1e-06
    ntax = 200, strht = 500000, srate = 1e-07
    ntax = 200, strht = 500000, srate = 1e-06
        Remove replicate 8
        Remove replicate 15
        Remove replicate 49
    """
    ntax = 200
    strhts = [10000000, 2000000, 500000]
    srates = [0.0000001, 0.000001]
    repls = range(1, 51)

    df = pandas.read_csv("true_vs_estimated_gene_trees.csv")

    for strht in strhts:
        for srate in srates:
            print("ntax = %d, strht = %g, srate = %g" % (ntax, strht, srate))
            for repl in repls:
                xdf = df[(df["NTAX"] == ntax) &
                         (df["STRHT"] == strht) &
                         (df["SRATE"] == srate) &
                         (df["REPL"] == repl)]

                nprune = 0
                for gene in range(1, 1001):
                    ydf = xdf[xdf["GENE"] == gene]

                    if ydf.shape[0] != 1:
                        sys.exit("ERROR!")

                    ni1 = float(ydf.I1.values[0])
                    ni2 = float(ydf.I2.values[0])

                    if (ni2 / ni1) < 0.50:
                        nprune += 1

                nkeep = 1000 - nprune
                if nkeep < 500:
                    print("    Remove replicate %d" % repl)

def find_bad_repls_non200tax():
    """
    ntax = 10, strht = 2e+06, srate = 1e-06
        Remove replicate 41
    ntax = 50, strht = 2e+06, srate = 1e-06
        Remove replicate 21
        Remove replicate 41
    ntax = 100, strht = 2e+06, srate = 1e-06
        Remove replicate 8
        Remove replicate 47
    ntax = 500, strht = 2e+06, srate = 1e-06
    ntax = 1000, strht = 2e+06, srate = 1e-06
    """
    ntaxs = [10, 50, 100, 500, 1000]
    strht = 2000000
    srate = 0.000001
    repls = range(1, 51)

    df = pandas.read_csv("true_vs_estimated_gene_trees.csv")

    for ntax in ntaxs:
        print("ntax = %d, strht = %g, srate = %g" % (ntax, strht, srate))
        for repl in repls:
            xdf = df[(df["NTAX"] == ntax) &
                     (df["STRHT"] == strht) &
                     (df["SRATE"] == srate) &
                     (df["REPL"] == repl)]

            nprune = 0
            for gene in range(1, 1001):
                ydf = xdf[xdf["GENE"] == gene]

                if ydf.shape[0] != 1:
                    sys.exit("ERROR!")

                ni1 = float(ydf.I1.values[0])
                ni2 = float(ydf.I2.values[0])

                if (ni2 / ni1) < 0.50:
                    nprune += 1

            nkeep = 1000 - nprune
            #print("  repl %d - %d" % (repl, nkeep))
            if nkeep < 500:
                print("    Remove replicate %d" % repl)

def report_stats_200tax():
    """
    10000000 1e-07 : GTEE = 0.283781 +/- 0.111238
    10000000 1e-06 : GTEE = 0.212423 +/- 0.124530
    2000000 1e-07 : GTEE = 0.334524 +/- 0.118986
    2000000 1e-06 : GTEE = 0.265774 +/- 0.122493
    500000 1e-07 : GTEE = 0.440191 +/- 0.140508
    500000 1e-06 : GTEE = 0.438317 +/- 0.121045

    10000000 1e-07 : ILS = 0.089425 +/- 0.010128
    10000000 1e-06 : ILS = 0.212668 +/- 0.024095
    2000000 1e-07 : ILS = 0.338618 +/- 0.016798
    2000000 1e-06 : ILS = 0.337251 +/- 0.018133
    500000 1e-07 : ILS = 0.683870 +/- 0.015363
    500000 1e-06 : ILS = 0.687386 +/- 0.016732

    10000000 1e-07 : AD = 0.305523 +/- 0.104339
    10000000 1e-06 : AD = 0.326050 +/- 0.093143
    2000000 1e-07 : AD = 0.466105 +/- 0.078986
    2000000 1e-06 : AD = 0.441597 +/- 0.073647
    500000 1e-07 : AD = 0.734598 +/- 0.027473
    500000 1e-06 : AD = 0.736384 +/- 0.030822
    """
    ntax = 200
    strhts = [10000000, 2000000, 500000]
    srates = [0.0000001, 0.000001]

    # GTEE
    df = pandas.read_csv("true_vs_estimated_gene_trees.csv")
    for strht in strhts:
        for srate in srates:
            rfs = []

            repls = list(range(1, 51))
            if (strht == 500000) and (srate == 0.000001):
                repls.remove(8)
                repls.remove(15)
                repls.remove(49)

            
            for repl in repls:
                xdf = df[(df["NTAX"] == ntax) &
                         (df["STRHT"] == strht) &
                         (df["SRATE"] == srate) &
                         (df["REPL"] == repl)]

                rfs.append(numpy.mean(xdf.RF.values))
            x = numpy.mean(rfs)
            y = numpy.std(rfs)
            print("%d %g : GTEE = %f +/- %f" % (strht, srate, x, y))


    # ILS
    df = pandas.read_csv("true_species_tree_vs_true_gene_trees.csv")
    for strht in strhts:
        for srate in srates:
            rfs = []

            repls = list(range(1, 51))
            if (strht == 500000) and (srate == 0.000001):
                repls.remove(8)
                repls.remove(15)
                repls.remove(49)
            
            for repl in repls:
                xdf = df[(df["NTAX"] == ntax) &
                         (df["STRHT"] == strht) &
                         (df["SRATE"] == srate) &
                         (df["REPL"] == repl)]

                rfs.append(numpy.mean(xdf.RF.values))

            x = numpy.mean(rfs)
            y = numpy.std(rfs)
            print("%d %g : ILS = %f +/- %f" % (strht, srate, x, y))

    # AD
    df = pandas.read_csv("true_species_tree_vs_estimated_gene_trees.csv")
    for strht in strhts:
        for srate in srates:
            rfs = []

            repls = list(range(1, 51))
            if (strht == 500000) and (srate == 0.000001):
                repls.remove(8)
                repls.remove(15)
                repls.remove(49)
            
            for repl in repls:
                xdf = df[(df["NTAX"] == ntax) &
                         (df["STRHT"] == strht) &
                         (df["SRATE"] == srate) &
                         (df["REPL"] == repl)]

                rfs.append(numpy.mean(xdf.RF.values))

            x = numpy.mean(rfs)
            y = numpy.std(rfs)
            print("%d %g : AD = %f +/- %f" % (strht, srate, x, y))


def report_stats_non200tax():
    """
    10 2000000 1e-06 : GTEE = 0.194153 +/- 0.086913
    50 2000000 1e-06 : GTEE = 0.256299 +/- 0.112229
    100 2000000 1e-06 : GTEE = 0.261846 +/- 0.090953
    500 2000000 1e-06 : GTEE = 0.283523 +/- 0.108514
    1000 2000000 1e-06 : GTEE = 0.297056 +/- 0.109003

    10 2000000 1e-06 : ILS = 0.171061 +/- 0.055135
    50 2000000 1e-06 : ILS = 0.306811 +/- 0.038339
    100 2000000 1e-06 : ILS = 0.332448 +/- 0.023935
    500 2000000 1e-06 : ILS = 0.342649 +/- 0.010347
    1000 2000000 1e-06 : ILS = 0.346056 +/- 0.008473

    10 2000000 1e-06 : AD = 0.281852 +/- 0.080144
    50 2000000 1e-06 : AD = 0.417533 +/- 0.076419
    100 2000000 1e-06 : AD = 0.435805 +/- 0.054340
    500 2000000 1e-06 : AD = 0.452495 +/- 0.065103
    1000 2000000 1e-06 : AD = 0.470366 +/- 0.075857
    """
    ntaxs = [10, 50, 100, 500, 1000]
    strht = 2000000
    srate = 0.000001
    repls = range(1, 51)

    # GTEE
    df = pandas.read_csv("true_vs_estimated_gene_trees.csv")
    for ntax in ntaxs:
        rfs = []

        repls = list(range(1, 51))
        if ntax == 10:
            repls.remove(41)   
        elif ntax == 50:
            repls.remove(21)
            repls.remove(41)
        elif ntax == 100:
            repls.remove(8)
            repls.remove(47)
        elif ntax == 1000:
            # Because ASTRAL-III failed to complete
            repls.remove(6)
            repls.remove(8)
            repls.remove(38)
            
        for repl in repls:
            xdf = df[(df["NTAX"] == ntax) &
                     (df["STRHT"] == strht) &
                     (df["SRATE"] == srate) &
                     (df["REPL"] == repl)]

            rfs.append(numpy.mean(xdf.RF.values))
        x = numpy.mean(rfs)
        y = numpy.std(rfs)
        print("%d %d %g : GTEE = %f +/- %f" % (ntax, strht, srate, x, y))


    # ILS
    df = pandas.read_csv("true_species_tree_vs_true_gene_trees.csv")
    for ntax in ntaxs:
        rfs = []

        repls = list(range(1, 51))
        if ntax == 10:
            repls.remove(41)   
        elif ntax == 50:
            repls.remove(21)
            repls.remove(41)
        elif ntax == 100:
            repls.remove(8)
            repls.remove(47)

        for repl in repls:
            xdf = df[(df["NTAX"] == ntax) &
                     (df["STRHT"] == strht) &
                     (df["SRATE"] == srate) &
                     (df["REPL"] == repl)]

            rfs.append(numpy.mean(xdf.RF.values))

        x = numpy.mean(rfs)
        y = numpy.std(rfs)
        print("%d %d %g : ILS = %f +/- %f" % (ntax, strht, srate, x, y))

    # AD
    df = pandas.read_csv("true_species_tree_vs_estimated_gene_trees.csv")
    for ntax in ntaxs:
        rfs = []

        repls = list(range(1, 51))
        if ntax == 10:
            repls.remove(41)   
        elif ntax == 50:
            repls.remove(21)
            repls.remove(41)
        elif ntax == 100:
            repls.remove(8)
            repls.remove(47)

        for repl in repls:
            xdf = df[(df["NTAX"] == ntax) &
                     (df["STRHT"] == strht) &
                     (df["SRATE"] == srate) &
                     (df["REPL"] == repl)]

            rfs.append(numpy.mean(xdf.RF.values))

        x = numpy.mean(rfs)
        y = numpy.std(rfs)
        print("%d %d %g : AD = %f +/- %f" % (ntax, strht, srate, x, y))


if __name__ == "__main__":
    #find_bad_repls_200tax()
    #find_bad_repls_non200tax()
    report_stats_200tax()
    report_stats_non200tax()
