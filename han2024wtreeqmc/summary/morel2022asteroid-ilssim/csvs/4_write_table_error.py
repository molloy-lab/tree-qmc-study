import pandas
import numpy
import sys

#sys.exit("DONE RUNNING")

def report_stats(df, mthds, ntax, ngen, nbps, blsc, psiz, miss, vary, emet, doit):
    sys.stdout.write("& %s" % vary)

    xdf = df[(df["NTAX"] == ntax) &
             (df["NGEN"] == ngen) &
             (df["NBPS"] == nbps) &
             (df["BLSC"] == blsc) &
             (df["PSIZ"] == psiz) &
             (df["MISS"] == miss)]

    keep = []
    minval = 1
    for ind, mthd in enumerate(mthds):
        if emet == "FNR":
            ers = xdf[(xdf["MTHD"] == mthd)].SEFNR.values
        elif emet == "FPR":
            ers = xdf[(xdf["MTHD"] == mthd)].SEFPR.values
        else:
            sys.exit("unknown error metric\n")

        if len(ers) != 50:
            print(ers)
            sys.exit("ERROR")

        x = round(round(numpy.mean(ers), 5), 4)
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

    if ((doit != "varypsiz") and \
        (ntax == 50) and \
        (ngen == 1000) and \
        (nbps == 100) and \
        (blsc == 1.0) and \
        (psiz == 50000000) and \
        (miss == 0.6)):
        sys.stdout.write(" (dup)")

    sys.stdout.write(" \\\\\n")


if __name__ == "__main__":
    
    emet = "FNR"
    emet = "FPR"

    sys.stdout.write("\\begin{table}[!h]\n")
    sys.stdout.write("\\caption[")
    if emet == "FNR":
        sys.stdout.write("Species tree error (FN) for Asteroid data]{\\textbf{Species tree (FN) error rate for Asteroid data.} ")
    else:
        sys.stdout.write("Species tree error (FP) for Asteroid data]{\\textbf{Species tree (FP) error rate for Asteroid data.} ")
    sys.stdout.write("Mean error rate is given across 50 replicates for each method.}\n")
    #sys.stdout.write("\\label{tab:}\n")
    sys.stdout.write("\\centering\n")
    sys.stdout.write("\\footnotesize\n")
    #sys.stdout.write("\\begin{tabular}{r c c c c c c c l}\n")
    sys.stdout.write("\\begin{tabular}{r c c c c c c c l}\n")
    sys.stdout.write("\\toprule \n")

    sys.stdout.write(" & & ASTRID & ASTER & Asteroid & TREE-QMC-n2 & n1 & n0 & n2-shared \\\\\n")
    sys.stdout.write("\\midrule\n")

    experiments = ["varypsiz", "varyntax", "varyngen",
                   "varynbps", "varyblsc", "varymiss"]

    mthds = ["wastrid_vanilla",
             "aster_v1.16.3.4",
             "asteroid",
             "wtreeqmc_wf_n2",
             "wtreeqmc_wf_n1",
             "wtreeqmc_wf_n0",
             "wtreeqmc_wf_n2_shared"] #, "wtreeqmc_wf_n1_shared"]

    # Sanity check by refining TREE-QMC trees
    #mthds = ["wastrid_vanilla",
    #         "aster_v1.16.3.4",
    #         "asteroid",
    #         "wtreeqmc_wf_n2_refined",
    #         "wtreeqmc_wf_n1_refined",
    #         "wtreeqmc_wf_n0_refined",
    #         "wtreeqmc_wf_n2_shared_refined"] #, "wtreeqmc_wf_n1_shared_refined"]

    namemap = {}
    namemap["asteroid"] = "Asteroid"
    namemap["aster_v1.16.3.4"] = "ASTER"
    namemap["wastrid_vanilla"] = "ASTRID"
    namemap["wtreeqmc_wf_n0"] = "TREE-QMC-n0"
    namemap["wtreeqmc_wf_n1"] = "TREE-QMC-n1"
    namemap["wtreeqmc_wf_n1_shared"] = "TREE-QMC-n1_shared"
    namemap["wtreeqmc_wf_n2"] = "TREE-QMC-n2"
    namemap["wtreeqmc_wf_n2_shared"] = "TREE-QMC-n2_shared"    

    for doit in experiments:
        sys.stdout.write("\\midrule\n")

        ntax = 50
        ngen = 1000
        nbps = 100
        blsc = 1.0
        psiz = 50000000
        miss = 0.6  # ms always same as mf

        df = pandas.read_csv("data-" + doit + "-error-and-timings.csv",
                             na_values='NA', keep_default_na=False)

        if doit == "varypsiz":
            #sys.stdout.write("\\multicolumn{5}{l}{\\textit{Varying population size}} \\\\[0.25em]\n")
            sys.stdout.write("\\multirow{ 5}{2cm}{Varying population size}\n")
            names = ["        10", 
                     "  50000000",
                     " 100000000",
                     " 500000000",
                     "1000000000"]
            psizs = [10, 50000000, 100000000, 500000000, 1000000000]
            for name, psiz in zip(names, psizs):
                report_stats(df, mthds, ntax, ngen, nbps, blsc, psiz, miss, name, emet, doit)
        elif doit == "varyntax":
            #sys.stdout.write("\\multicolumn{5}{l}{\\textit{Varying number of taxa}} \\\\[0.25em]\n")
            sys.stdout.write("\\multirow{ 6}{2cm}{Varying number of taxa}\n")
            names = ["  25", "  75", "  50", " 100", " 125", " 150"]
            ntaxs = [25, 75, 50, 100, 125, 150]  # dup 50
            for name, ntax in zip(names, ntaxs):
                report_stats(df, mthds, ntax, ngen, nbps, blsc, psiz, miss, name, emet, doit)
        elif doit == "varyngen":
            #sys.stdout.write("\\multicolumn{5}{l}{\\textit{Varying number of genes}} \\\\[0.25em]\n")
            sys.stdout.write("\\multirow{ 4}{2cm}{Varying number of genes}\n")
            names = [" 250", " 500", "1000", "2000"]
            ngens = [250, 500, 1000, 2000]  # dup 1000
            for name, ngen in zip(names, ngens):
                report_stats(df, mthds, ntax, ngen, nbps, blsc, psiz, miss, name, emet, doit)
        elif doit == "varynbps":
            #sys.stdout.write("\\multicolumn{5}{l}{\\textit{Varying sequence length}} \\\\[0.25em]\n")
            sys.stdout.write("\\multirow{ 4}{2cm}{Varying sequence length}\n")
            names = ["  50", " 100", " 200", " 500"]
            nbpss = [50, 100, 200, 500]  # dup 100
            for name, nbps in zip(names, nbpss):
                report_stats(df, mthds, ntax, ngen, nbps, blsc, psiz, miss, name, emet, doit)
        elif doit == "varyblsc":
            #sys.stdout.write("\\multicolumn{5}{l}{\\textit{Varying branch length scaler}} \\\\[0.25em]\n")
            sys.stdout.write("\\multirow{ 6}{2cm}{Varying branch length scaler}\n")
            names = ["  0.05", "  0.10", "  1.00", " 10.00", "100.00", "200.00"]
            blscs = [0.05, 0.1, 1, 10, 100, 200]  # dup 1.0
            for name, blsc in zip(names, blscs):
                report_stats(df, mthds, ntax, ngen, nbps, blsc, psiz, miss, name, emet, doit)
        elif doit == "varymiss":
            #sys.stdout.write("\\multicolumn{5}{l}{\\textit{Varying missingness parameter}} \\\\[0.25em]\n")
            sys.stdout.write("\\multirow{ 6}{2cm}{Varying missingness parameter}\n")
            names = ["0.50", "0.55", "0.60", "0.65", "0.70", "0.75"]
            misss = [0.5, 0.55, 0.6, 0.65, 0.7, 0.75]  # dup 0.6
            for name, miss in zip(names, misss):
                report_stats(df, mthds, ntax, ngen, nbps, blsc, psiz, miss, name, emet, doit)
        else:
            sys.exit()

    sys.stdout.write("\\bottomrule\n")
    sys.stdout.write("\\end{tabular}\n")
    sys.stdout.write("\\end{table}\n")


"""
\begin{table}[!h]
\caption[Species tree error (FP) for Asteroid data]{\textbf{Species tree (FP) error rate for Asteroid data.} Mean error rate is given across 50 replicates for each method.}
\centering
\footnotesize
\begin{tabular}{r c c c c c c c l}
\toprule 
 & & ASTRID & ASTER & Asteroid & TREE-QMC-n2 & n1 & n0 & n2-shared \\
\midrule
\midrule
\multirow{ 5}{2cm}{Varying population size}
&         10 & 0.21200 & 0.17660 & 0.14430 & \textbf{0.14370} & 0.14630 & 0.17560 & 0.21490 \\
&   50000000 & 0.24600 & 0.20380 & 0.18130 & \textbf{0.18120} & 0.18190 & 0.20710 & 0.24580 \\
&  100000000 & 0.29660 & 0.26170 & 0.23150 & \textbf{0.22390} & 0.22740 & 0.25870 & 0.32100 \\
&  500000000 & 0.47530 & 0.46210 & 0.43450 & \textbf{0.41560} & 0.42080 & 0.45200 & 0.51840 \\
& 1000000000 & 0.58470 & 0.58810 & \textbf{0.55280} & 0.55330 & 0.55430 & 0.60310 & 0.66180 \\
\midrule
\multirow{ 6}{2cm}{Varying number of taxa}
&   25 & 0.26910 & 0.27540 & 0.22540 & 0.22280 & \textbf{0.22160} & 0.26480 & 0.27030 \\
&   75 & 0.24390 & 0.20890 & 0.18310 & \textbf{0.18070} & 0.18420 & 0.20710 & 0.26900 \\
&   50 & 0.24600 & 0.20380 & 0.18130 & \textbf{0.18120} & 0.18190 & 0.20710 & 0.24580 (dup) \\
&  100 & 0.25070 & 0.21710 & 0.19650 & \textbf{0.18560} & 0.18980 & 0.22000 & 0.28340 \\
&  125 & 0.24670 & 0.21120 & 0.18480 & \textbf{0.17650} & 0.18010 & 0.20990 & 0.27040 \\
&  150 & 0.25580 & 0.21690 & 0.18680 & \textbf{0.18160} & 0.18600 & 0.21720 & 0.28880 \\
\midrule
\multirow{ 4}{2cm}{Varying number of genes}
&  250 & 0.47080 & 0.43140 & \textbf{0.35980} & 0.36450 & 0.36500 & 0.42440 & 0.46060 \\
&  500 & 0.34820 & 0.31340 & 0.26010 & \textbf{0.25260} & 0.25960 & 0.29610 & 0.33560 \\
& 1000 & 0.24600 & 0.20380 & 0.18130 & \textbf{0.18120} & 0.18190 & 0.20710 & 0.24580 (dup) \\
& 2000 & 0.18170 & 0.14940 & 0.13110 & \textbf{0.12050} & 0.12630 & 0.14860 & 0.18800 \\
\midrule
\multirow{ 4}{2cm}{Varying sequence length}
&   50 & 0.29180 & 0.29470 & \textbf{0.24580} & 0.25020 & 0.25440 & 0.29120 & 0.34300 \\
&  100 & 0.24600 & 0.20380 & 0.18130 & \textbf{0.18120} & 0.18190 & 0.20710 & 0.24580 (dup) \\
&  200 & 0.23060 & 0.17110 & 0.16430 & \textbf{0.15320} & 0.15450 & 0.17300 & 0.20780 \\
&  500 & 0.19740 & 0.14640 & 0.13230 & 0.11520 & \textbf{0.11170} & 0.13870 & 0.15590 \\
\midrule
\multirow{ 6}{2cm}{Varying branch length scaler}
&   0.05 & 0.51020 & 0.48550 & 0.45060 & 0.45540 & \textbf{0.44880} & 0.48060 & 0.58390 \\
&   0.10 & 0.37830 & 0.34000 & \textbf{0.29280} & 0.30920 & 0.30210 & 0.35010 & 0.42360 \\
&   1.00 & 0.24600 & 0.20380 & 0.18130 & \textbf{0.18120} & 0.18190 & 0.20710 & 0.24580 (dup) \\
&  10.00 & 0.26130 & 0.24510 & 0.20720 & \textbf{0.19000} & 0.20270 & 0.23790 & 0.26610 \\
& 100.00 & 0.42680 & 0.42550 & \textbf{0.35790} & 0.36260 & 0.37720 & 0.42630 & 0.51750 \\
& 200.00 & 0.50080 & 0.49530 & \textbf{0.45620} & 0.48220 & 0.48530 & 0.52540 & 0.61540 \\
\midrule
\multirow{ 6}{2cm}{Varying missingness parameter}
& 0.50 & 0.13360 & 0.11790 & \textbf{0.09870} & 0.10280 & 0.10490 & 0.12080 & 0.14070 \\
& 0.55 & 0.17320 & 0.17190 & 0.14510 & \textbf{0.14020} & 0.14300 & 0.17010 & 0.20660 \\
& 0.60 & 0.24600 & 0.20380 & 0.18130 & \textbf{0.18120} & 0.18190 & 0.20710 & 0.24580 (dup) \\
& 0.65 & 0.38490 & 0.32350 & 0.28050 & \textbf{0.24730} & 0.25690 & 0.30250 & 0.35640 \\
& 0.70 & 0.60670 & 0.49840 & 0.43100 & \textbf{0.41460} & 0.41930 & 0.45520 & 0.53320 \\
& 0.75 & 0.82050 & 0.73420 & 0.62240 & 0.61020 & \textbf{0.60340} & 0.65720 & 0.72410 \\
\bottomrule
\end{tabular}
\end{table}

\begin{table}[!h]
\caption[Species tree error (FN) for Asteroid data]{\textbf{Species tree (FN) error rate for Asteroid data.} Mean error rate is given across 50 replicates for each method.}
\centering
\footnotesize
\begin{tabular}{r c c c c c c c l}
\toprule 
 & & ASTRID & ASTER & Asteroid & TREE-QMC-n2 & n1 & n0 & n2-shared \\
\midrule
\midrule
\multirow{ 5}{2cm}{Varying population size}
&         10 & 0.21200 & 0.17660 & \textbf{0.14430} & 0.15750 & 0.16090 & 0.18860 & 0.22640 \\
&   50000000 & 0.24600 & 0.20380 & \textbf{0.18130} & 0.19620 & 0.19620 & 0.22040 & 0.25830 \\
&  100000000 & 0.29660 & 0.26170 & \textbf{0.23150} & 0.23960 & 0.24340 & 0.27400 & 0.33400 \\
&  500000000 & 0.47530 & 0.46210 & 0.43450 & \textbf{0.42430} & 0.42980 & 0.46080 & 0.52430 \\
& 1000000000 & 0.58470 & 0.58810 & \textbf{0.55280} & 0.56130 & 0.56170 & 0.61020 & 0.66890 \\
\midrule
\multirow{ 6}{2cm}{Varying number of taxa}
&   25 & 0.26910 & 0.27540 & \textbf{0.22540} & 0.24180 & 0.24000 & 0.28090 & 0.28820 \\
&   75 & 0.24390 & 0.20890 & \textbf{0.18310} & 0.19560 & 0.19970 & 0.22190 & 0.28420 \\
&   50 & 0.24600 & 0.20380 & \textbf{0.18130} & 0.19620 & 0.19620 & 0.22040 & 0.25830 (dup) \\
&  100 & 0.25070 & 0.21710 & \textbf{0.19650} & 0.20310 & 0.20640 & 0.23570 & 0.29770 \\
&  125 & 0.24670 & 0.21120 & \textbf{0.18480} & 0.19430 & 0.19900 & 0.22590 & 0.28560 \\
&  150 & 0.25580 & 0.21690 & \textbf{0.18680} & 0.19820 & 0.20250 & 0.23100 & 0.30290 \\
\midrule
\multirow{ 4}{2cm}{Varying number of genes}
&  250 & 0.47080 & 0.43140 & \textbf{0.35980} & 0.41230 & 0.41310 & 0.46560 & 0.50060 \\
&  500 & 0.34820 & 0.31340 & \textbf{0.26010} & 0.28310 & 0.29040 & 0.32360 & 0.36020 \\
& 1000 & 0.24600 & 0.20380 & \textbf{0.18130} & 0.19620 & 0.19620 & 0.22040 & 0.25830 (dup) \\
& 2000 & 0.18170 & 0.14940 & 0.13110 & \textbf{0.13060} & 0.13570 & 0.15830 & 0.19700 \\
\midrule
\multirow{ 4}{2cm}{Varying sequence length}
&   50 & 0.29180 & 0.29470 & \textbf{0.24580} & 0.26320 & 0.26790 & 0.30100 & 0.35390 \\
&  100 & 0.24600 & 0.20380 & \textbf{0.18130} & 0.19620 & 0.19620 & 0.22040 & 0.25830 (dup) \\
&  200 & 0.23060 & 0.17110 & \textbf{0.16430} & 0.17530 & 0.17570 & 0.19280 & 0.22720 \\
&  500 & 0.19740 & 0.14640 & \textbf{0.13230} & 0.13740 & 0.13400 & 0.15960 & 0.17570 \\
\midrule
\multirow{ 6}{2cm}{Varying branch length scaler}
&   0.05 & 0.51020 & 0.48550 & \textbf{0.45060} & 0.46850 & 0.46300 & 0.48980 & 0.59190 \\
&   0.10 & 0.37830 & 0.34000 & \textbf{0.29280} & 0.32260 & 0.31570 & 0.36080 & 0.43360 \\
&   1.00 & 0.24600 & 0.20380 & \textbf{0.18130} & 0.19620 & 0.19620 & 0.22040 & 0.25830 (dup) \\
&  10.00 & 0.26130 & 0.24510 & \textbf{0.20720} & 0.21150 & 0.22260 & 0.25570 & 0.28380 \\
& 100.00 & 0.42680 & 0.42550 & \textbf{0.35790} & 0.37230 & 0.38770 & 0.43450 & 0.52770 \\
& 200.00 & 0.50080 & 0.49530 & \textbf{0.45620} & 0.49320 & 0.49490 & 0.53280 & 0.62340 \\
\midrule
\multirow{ 6}{2cm}{Varying missingness parameter}
& 0.50 & 0.13360 & 0.11790 & \textbf{0.09870} & 0.10470 & 0.10680 & 0.12300 & 0.14260 \\
& 0.55 & 0.17320 & 0.17190 & \textbf{0.14510} & 0.14720 & 0.14940 & 0.17450 & 0.21280 \\
& 0.60 & 0.24600 & 0.20380 & \textbf{0.18130} & 0.19620 & 0.19620 & 0.22040 & 0.25830 (dup) \\
& 0.65 & 0.38490 & 0.32350 & \textbf{0.28050} & 0.29030 & 0.29970 & 0.33970 & 0.38820 \\
& 0.70 & 0.60670 & 0.49840 & \textbf{0.43100} & 0.48830 & 0.49130 & 0.51890 & 0.58580 \\
& 0.75 & 0.82050 & 0.73420 & \textbf{0.62240} & 0.69620 & 0.69020 & 0.72750 & 0.77880 \\
\bottomrule
\end{tabular}
\end{table}
"""