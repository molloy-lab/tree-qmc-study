import pandas
import numpy
import sys

sys.exit("DONE RUNNING")

def report_stats(df):
    repls = range(3000, 3050)

    keep = []

    repls = list(set(df.SEED.values))
    for repl in repls:
        xdf = df[df["SEED"] == repl]
        fnrs = xdf.FN.values / xdf.I1.values
        keep.append(numpy.mean(fnrs))

    x = round(round(numpy.mean(keep), 5), 4)
    y = round(round(numpy.std(keep), 5), 4)
    sys.stdout.write(" & $%1.4f \\pm %1.4f$" % (x, y))
    #sys.stdout.write(" & %1.4f" % (x))

def report_stats_miss(df, ntax):
    repls = range(3000, 3050)

    keep = []

    repls = list(set(df.SEED.values))
    for repl in repls:
        xdf = df[df["SEED"] == repl]
        miss = 1.0 - (xdf.NL.values / ntax)
        keep.append(numpy.mean(miss))

    x = round(round(numpy.mean(keep), 5), 4)
    y = round(round(numpy.std(keep), 5), 4)
    sys.stdout.write(" & $%1.4f \\pm %1.4f$" % (x, y))
    #sys.stdout.write(" & %1.4f" % (x))

def run_stats_helper(df_ils, df_gtee, df_ad, ntax, ngen, nbps, blsc, psiz, miss, vary, doit):
    sys.stdout.write("& %s" % vary)

    xdf_ils = df_ils[(df_ils["S"] == ntax) & \
                     (df_ils["F"] == ngen) & \
                     (df_ils["SITES"] == nbps) &
                     (df_ils["BL"] == blsc) & \
                     (df_ils["POP"] == psiz) & \
                     (df_ils["MS"] == miss)]
    report_stats(xdf_ils)

    xdf_gtee = df_gtee[(df_gtee["S"] == ntax) & \
                       (df_gtee["F"] == ngen) & \
                       (df_gtee["SITES"] == nbps) &
                       (df_gtee["BL"] == blsc) & \
                       (df_gtee["POP"] == psiz) & \
                       (df_gtee["MS"] == miss)]
    report_stats(xdf_gtee)

    xdf_ad = df_ad[(df_ad["S"] == ntax) & \
                   (df_ad["F"] == ngen) & \
                   (df_ad["SITES"] == nbps) &
                   (df_ad["BL"] == blsc) & \
                   (df_ad["POP"] == psiz) & \
                   (df_ad["MS"] == miss)]
    report_stats(xdf_ad)

    report_stats_miss(xdf_ad, int(ntax.replace('s', '')))

    if ((doit != "varypsiz") and \
        (ntax == "s50") and \
        (ngen == "f1000") and \
        (nbps == "sites100") and \
        (blsc == "bl1.0") and \
        (psiz == "pop50000000") and \
        (miss == "ms0.6")):
        sys.stdout.write(" (dup)")

    sys.stdout.write(" \\\\\n")

if __name__ == "__main__":
    df_ils = pandas.read_csv("all_true_species_tree_vs_true_gene_trees.csv")
    df_gtee = pandas.read_csv("all_true_vs_estimated_gene_trees.csv")
    df_ad = pandas.read_csv("all_true_species_tree_vs_estimated_gene_trees.csv")

    sys.stdout.write("\\begin{table}[!h]\n")
    sys.stdout.write("\\caption[Properties of Asteroid simulated data]")
    sys.stdout.write("{\\textbf{Properties of Asteroid simulated data from \\cite{morel2022asteroid}.} Data were downloaded from ")
    sys.stdout.write("\\url{https://cme.h-its.org/exelixis/material/asteroid_data.tar.gz} in December 2022. ")
    sys.stdout.write("For each replicate, ILS is the fraction of branches in true species tree missing from the true gene tree, averaged across all gene trees. ")
    sys.stdout.write("GTEE is the fraction of branches in true gene tree that are missing from the estimated gene tree, averaged across all gene trees. ")
    sys.stdout.write("AD is the fraction of branches in the true species tree that are missing from the estimated gene tree, averaged across all gene trees. ") 
    sys.stdout.write("The values in the table are the average ($\\pm$ standard deviation) across all replicates.}\n")
    sys.stdout.write("\\label{tab:morel}\n")
    sys.stdout.write("\\centering\n")
    sys.stdout.write("\\footnotesize\n")
    sys.stdout.write("\\begin{tabular}{l r c c c l}\n") # only need if using multirow 
    #sys.stdout.write("\\begin{tabular}{r c c c l}\n")
    sys.stdout.write("\\toprule \n")
    sys.stdout.write("\\multirow{ 1}{2cm}{}\n")  # only need if using multirow
    sys.stdout.write(" & & ILS & GTEE & AD & Proportion missing taxa \\\\\n")

    experiments = ["varypsiz", "varyntax", "varyngen",
                   "varynbps", "varyblsc", "varymiss"]

    for doit in experiments:
        sys.stdout.write("\\midrule\n")

        ntax = "s50"
        ngen = "f1000"
        nbps = "sites100"
        blsc = "bl1.0"
        psiz = "pop50000000"
        miss = "ms0.6"  # always same as mf

        if doit == "varypsiz":
            #sys.stdout.write("\\multicolumn{6}{l}{\\textit{Varying population size}} \\\\[0.25em]\n")
            sys.stdout.write("\\multirow{ 5}{2cm}{Varying population size}\n")
            names = ["        10", 
                     "  50000000",
                     " 100000000",
                     " 500000000",
                     "1000000000"]
            psizs = ["pop10", "pop50000000", "pop100000000", "pop500000000", "pop1000000000"]
            for name, psiz in zip(names, psizs):
                run_stats_helper(df_ils, df_gtee, df_ad,
                                 ntax, ngen, nbps, blsc, psiz, miss,
                                 name, doit)
        elif doit == "varyntax":
            #sys.stdout.write("\\multicolumn{6}{l}{\\textit{Varying number of taxa}} \\\\[0.25em]\n")
            sys.stdout.write("\\multirow{ 6}{2cm}{Varying number of taxa}\n")
            names = ["  25", "  75", "  50", " 100", " 125", " 150"]
            ntaxs = ["s25", "s75", "s50", "s100", "s125", "s150"]
            for name, ntax in zip(names, ntaxs):
                run_stats_helper(df_ils, df_gtee, df_ad,
                                 ntax, ngen, nbps, blsc, psiz, miss,
                                 name, doit)
        elif doit == "varyngen":
            #sys.stdout.write("\\multicolumn{6}{l}{\\textit{Varying number of genes}} \\\\[0.25em]\n")
            sys.stdout.write("\\multirow{ 4}{2cm}{Varying number of genes}\n")
            names = [" 250", " 500", "1000", "2000"]
            ngens = ["f250", "f500", "f1000", "f2000"]
            for name, ngen in zip(names, ngens):
                run_stats_helper(df_ils, df_gtee, df_ad,
                                 ntax, ngen, nbps, blsc, psiz, miss,
                                 name, doit)
        elif doit == "varynbps":
            #sys.stdout.write("\\multicolumn{6}{l}{\\textit{Varying sequence length}} \\\\[0.25em]\n")
            sys.stdout.write("\\multirow{ 4}{2cm}{Varying sequence length}\n")
            names = ["  50", " 100", " 200", " 500"]
            nbpss = ["sites50", "sites100", "sites200", "sites500"]
            for name, nbps in zip(names, nbpss):
                run_stats_helper(df_ils, df_gtee, df_ad,
                                 ntax, ngen, nbps, blsc, psiz, miss,
                                 name, doit)
        elif doit == "varyblsc":
            #sys.stdout.write("\\multicolumn{6}{l}{\\textit{Varying branch length scaler}} \\\\[0.25em]\n")
            sys.stdout.write("\\multirow{ 6}{2cm}{Varying branch length scaler}\n")
            names = ["  0.05", "  0.10", "  1.00", " 10.00", "100.00", "200.00"]
            blscs = ["bl0.05", "bl0.1", "bl1.0", "bl10.0", "bl100.0", "bl200.0"]
            for name, blsc in zip(names, blscs):
                run_stats_helper(df_ils, df_gtee, df_ad,
                                 ntax, ngen, nbps, blsc, psiz, miss,
                                 name, doit)
        elif doit == "varymiss":
            #sys.stdout.write("\\multicolumn{6}{l}{\\textit{Varying missingness parameter}} \\\\[0.25em]\n")
            sys.stdout.write("\\multirow{ 6}{2cm}{Varying missingness parameter}\n")
            names = ["0.50", "0.55", "0.60", "0.65", "0.70", "0.75"]
            misss = ["ms0.5", "ms0.55", "ms0.6", "ms0.65", "ms0.7", "ms0.75"]
            for name, miss in zip(names, misss):
                run_stats_helper(df_ils, df_gtee, df_ad,
                                 ntax, ngen, nbps, blsc, psiz, miss,
                                 name, doit)
        else:
            sys.exit()

    sys.stdout.write("\\bottomrule\n")
    sys.stdout.write("\\end{tabular}\n")
    sys.stdout.write("\\end{table}\n")


"""
\begin{table}[!h]
\caption[Properties of Asteroid simulated data]{\textbf{Properties of Asteroid simulated data from \cite{morel2022asteroid}.} Data were downloaded from \url{https://cme.h-its.org/exelixis/material/asteroid_data.tar.gz} in December 2022. For each replicate, ILS is the fraction of branches in true species tree missing from the true gene tree, averaged across all gene trees. GTEE is the fraction of branches in true gene tree that are missing from the estimated gene tree, averaged across all gene trees. AD is the fraction of branches in the true species tree that are missing from the estimated gene tree, averaged across all gene trees. The values in the table are the average ($\pm$ standard deviation) across all replicates.}
\label{tab:morel}
\centering
\footnotesize
\begin{tabular}{l r c c c l}
\toprule 
\multirow{ 1}{2cm}{}
 & & ILS & GTEE & AD & Proportion missing taxa \\
\midrule
\multirow{ 5}{2cm}{Varying population size}
&         10 & $0.0000 \pm 0.0000$ & $0.3419 \pm 0.0453$ & $0.3419 \pm 0.0453$ & $0.8115 \pm 0.0107$ \\
&   50000000 & $0.0610 \pm 0.0297$ & $0.3530 \pm 0.0483$ & $0.3616 \pm 0.0496$ & $0.8106 \pm 0.0118$ \\
&  100000000 & $0.1119 \pm 0.0359$ & $0.3598 \pm 0.0440$ & $0.3836 \pm 0.0478$ & $0.8109 \pm 0.0115$ \\
&  500000000 & $0.3854 \pm 0.0722$ & $0.3904 \pm 0.0413$ & $0.5428 \pm 0.0582$ & $0.8131 \pm 0.0099$ \\
& 1000000000 & $0.5576 \pm 0.0626$ & $0.4071 \pm 0.0395$ & $0.6579 \pm 0.0501$ & $0.8118 \pm 0.0089$ \\
\midrule
\multirow{ 6}{2cm}{Varying number of taxa}
&   25 & $0.0492 \pm 0.0314$ & $0.3008 \pm 0.0456$ & $0.3085 \pm 0.0485$ & $0.7607 \pm 0.0104$ \\
&   75 & $0.0626 \pm 0.0225$ & $0.3872 \pm 0.0382$ & $0.3960 \pm 0.0388$ & $0.8255 \pm 0.0088$ \\
&   50 & $0.0610 \pm 0.0297$ & $0.3530 \pm 0.0483$ & $0.3616 \pm 0.0496$ & $0.8106 \pm 0.0118$ (dup row) \\
&  100 & $0.0733 \pm 0.0270$ & $0.4125 \pm 0.0442$ & $0.4220 \pm 0.0454$ & $0.8320 \pm 0.0076$ \\
&  125 & $0.0835 \pm 0.0274$ & $0.4384 \pm 0.0397$ & $0.4487 \pm 0.0417$ & $0.8340 \pm 0.0073$ \\
&  150 & $0.0841 \pm 0.0284$ & $0.4524 \pm 0.0362$ & $0.4625 \pm 0.0371$ & $0.8370 \pm 0.0073$ \\
\midrule
\multirow{ 4}{2cm}{Varying number of genes}
&  250 & $0.0602 \pm 0.0256$ & $0.3504 \pm 0.0487$ & $0.3584 \pm 0.0500$ & $0.8127 \pm 0.0100$ \\
&  500 & $0.0592 \pm 0.0264$ & $0.3527 \pm 0.0435$ & $0.3603 \pm 0.0448$ & $0.8110 \pm 0.0094$ \\
& 1000 & $0.0610 \pm 0.0297$ & $0.3530 \pm 0.0483$ & $0.3616 \pm 0.0496$ & $0.8106 \pm 0.0118$ (dup row) \\
& 2000 & $0.0618 \pm 0.0256$ & $0.3496 \pm 0.0464$ & $0.3582 \pm 0.0478$ & $0.8116 \pm 0.0101$ \\
\midrule
\multirow{ 4}{2cm}{Varying sequence length}
&   50 & $0.0618 \pm 0.0264$ & $0.4578 \pm 0.0452$ & $0.4628 \pm 0.0455$ & $0.8113 \pm 0.0093$ \\
&  100 & $0.0610 \pm 0.0297$ & $0.3530 \pm 0.0483$ & $0.3616 \pm 0.0496$ & $0.8106 \pm 0.0118$ (dup row) \\
&  200 & $0.0586 \pm 0.0282$ & $0.2514 \pm 0.0454$ & $0.2634 \pm 0.0480$ & $0.8098 \pm 0.0100$ \\
&  500 & $0.0614 \pm 0.0273$ & $0.1571 \pm 0.0374$ & $0.1779 \pm 0.0423$ & $0.8135 \pm 0.0100$ \\
\midrule
\multirow{ 6}{2cm}{Varying branch length scaler}
&   0.05 & $0.0600 \pm 0.0274$ & $0.5530 \pm 0.0610$ & $0.5559 \pm 0.0612$ & $0.8135 \pm 0.0093$ \\
&   0.10 & $0.0608 \pm 0.0282$ & $0.4686 \pm 0.0555$ & $0.4722 \pm 0.0550$ & $0.8118 \pm 0.0090$ \\
&   1.00 & $0.0610 \pm 0.0297$ & $0.3530 \pm 0.0483$ & $0.3616 \pm 0.0496$ & $0.8106 \pm 0.0118$ (dup row) \\
&  10.00 & $0.0600 \pm 0.0262$ & $0.4280 \pm 0.0412$ & $0.4347 \pm 0.0414$ & $0.8125 \pm 0.0097$ \\
& 100.00 & $0.0599 \pm 0.0270$ & $0.7035 \pm 0.0382$ & $0.7050 \pm 0.0378$ & $0.8126 \pm 0.0102$ \\
& 200.00 & $0.0604 \pm 0.0258$ & $0.7664 \pm 0.0322$ & $0.7674 \pm 0.0319$ & $0.8109 \pm 0.0098$ \\
\midrule
\multirow{ 6}{2cm}{Varying missingness parameter}
& 0.50 & $0.0767 \pm 0.0281$ & $0.3818 \pm 0.0430$ & $0.3936 \pm 0.0441$ & $0.7373 \pm 0.0111$ \\
& 0.55 & $0.0659 \pm 0.0265$ & $0.3644 \pm 0.0467$ & $0.3744 \pm 0.0481$ & $0.7792 \pm 0.0124$ \\
& 0.60 & $0.0610 \pm 0.0297$ & $0.3530 \pm 0.0483$ & $0.3616 \pm 0.0496$ & $0.8106 \pm 0.0118$ (dup row) \\
& 0.65 & $0.0548 \pm 0.0250$ & $0.3340 \pm 0.0449$ & $0.3418 \pm 0.0459$ & $0.8378 \pm 0.0088$ \\
& 0.70 & $0.0471 \pm 0.0250$ & $0.3200 \pm 0.0462$ & $0.3263 \pm 0.0468$ & $0.8617 \pm 0.0093$ \\
& 0.75 & $0.0447 \pm 0.0267$ & $0.3023 \pm 0.0538$ & $0.3066 \pm 0.0539$ & $0.8789 \pm 0.0048$ \\
\bottomrule
\end{tabular}
\end{table}
"""
