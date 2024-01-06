import pandas
import numpy
import sys
import warnings

def report_stats(df, mthds, ngen, ntax, hght, rate, supp):
    sys.stdout.write("%d & %s & %s & %d" % (ntax, hght, rate, ngen))

    nrepl = len(df[(df["MTHD"] == "CA-ML")].REPL.values)

    keep = []
    minval = 1
    for ind, mthd in enumerate(mthds):
        if (mthd == "CA-ML") or \
           (mthd == "TQMC-n2-origstudy") or \
           (mthd == "TQMC-n2"):
            xdf = df[(df["MTHD"] == mthd)]
        else:
            xdf = df[(df["MTHD"] == mthd) & (df["SUPP"] == supp)]

        vals = xdf.SEFNR.values

        if len(vals) != nrepl:
            sys.exit("\nERROR - inconsistent number of replicates")
        if len(vals) > 50:
            sys.exit("\nERROR - too many replicates")

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            fnravg = numpy.nanmean(vals)

        x = round(round(fnravg, 5), 4)
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

    sys.stdout.write(" \\\\\n")

if __name__ == "__main__":
    mthds = ["CA-ML",
             "ASTER-wh",
             "TQMC-wh_n2",
             "TQMC-n2",
             "TQMC-n2-origstudy",
             "ASTRID-ws"]

    supps = ["abayes", "sh"]
    ngens = [50, 200, 1000]
    ntaxs = [10, 50, 100, 200, 500, 1000]

    sys.stdout.write("\\begin{table}[!h]\n")
    sys.stdout.write("\\caption[")
    sys.stdout.write("RF error rate for ASTRAL-II simulated data]{\\textbf{RF error rate for ASTRAL-II simulated data.} ")
    sys.stdout.write("Mean error rate is given across replicates for each method. ")
    #sys.stdout.write("All weighted methods were given gene trees with \\textbf{%s} support. " % supp)
    sys.stdout.write("The TREE-QMC-n2 from the original study cannot handle polytomies in gene trees and thus refined polytomies due to identical sequences randomly. ")
    sys.stdout.write("The TREE-QMC-n2 from this study handles polytomies by refining them and giving the produced branches support 0 (all other branches have support 1). ")
    sys.stdout.write("Note: the original study did not run methods on data sets with 50 or 200 genes. }\n")
    #sys.stdout.write("\\label{tab:}\n")
    sys.stdout.write("\\centering\n")
    sys.stdout.write("\\scriptsize\n")
    sys.stdout.write("\\begin{tabular}{c c c c c c c c c c c c}\n")
    sys.stdout.write("\\toprule \n")

    sys.stdout.write(" \\# of & ILS & speciation & \\# of & CA-ML & ASTER & TREE-QMC & TREE-QMC & TREE-QMC & ASTRID \\\\\n")
    sys.stdout.write(" taxa & level & location & genes &  & wh & wh-n2 & n2 & n2 (original study) & ws \\\\\n")

    for supp in supps:
        sys.stdout.write("\\midrule\n")
        sys.stdout.write("\\multicolumn{10}{c}{\\textbf{%s support for weighted methods}} \\\\\n" % supp)
        sys.stdout.write("\\midrule\n")
        for ntax in ntaxs:
            if ntax == 200:
                df = pandas.read_csv("data-varyils-error.csv", keep_default_na=True)
                hghts = ["low", "medium", "high"]
                rates = ["deep", "shallow"]
            else:
                df = pandas.read_csv("data-varyntax-error.csv", keep_default_na=True)
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

