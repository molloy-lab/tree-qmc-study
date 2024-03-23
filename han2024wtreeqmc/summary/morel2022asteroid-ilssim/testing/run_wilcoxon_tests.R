#!/usr/bin/env Rscript

library(exactRankTests)
library(coin)

run_tests_helper <- function(df, mthd, ntax, ngen, nbps, blsc, psiz, miss, name, emet, doit, threshold_bonferoni) {
    #print(mthd)
    #print(ntax)
    #print(ngen)
    #print(nbps)
    #print(blsc)
    #print(psiz)
    #print(miss)
    #print(name)
    #print(emet)
    #print(doit)

    xdf <- df[df$NTAX == ntax, ]
    xdf <- xdf[xdf$NGEN == ngen, ]
    xdf <- xdf[xdf$NBPS == nbps, ]
    xdf <- xdf[xdf$BLSC == blsc, ]
    xdf <- xdf[xdf$PSIZ == psiz, ]
    xdf <- xdf[xdf$MISS == miss, ]

    # Find ties
    if (emet == "FNR") { 
        mthd1 <- xdf$TQMCn2xFNR
        if (mthd == "ASTER") { mthd2 <- xdf$ASTERxFNR } 
        else if (mthd == "ASTRID") { mthd2 <- xdf$WASTRIDxFNR } 
        else if (mthd == "ASTEROID") { mthd2 <- xdf$ASTEROIDxFNR }
        else if (mthd == "TREE-QMC-n2-shared") { mthd2 <- xdf$TQMCn2sharedxFNR }
        else if (mthd == "TREE-QMC-n1") { mthd2 <- xdf$TQMCn1xFNR }
        else if (mthd == "TREE-QMC-n0") { mthd2 <- xdf$TQMCn0xFNR }
        else {
            print("ERROR")
            exit()
        }
    } else if (emet == "FPR") {
        mthd1 <- xdf$TQMCn2xFPR
        if (mthd == "ASTER") { mthd2 <- xdf$ASTERxFPR }
        else if (mthd == "ASTRID") { mthd2 <- xdf$WASTRIDxFPR }
        else if (mthd == "ASTEROID") { mthd2 <- xdf$ASTEROIDxFPR }
        else if (mthd == "TREE-QMC-n2-shared") { mthd2 <- xdf$TQMCn2sharedxFPR }
        else if (mthd == "TREE-QMC-n1") { mthd2 <- xdf$TQMCn1xFPR } 
        else if (mthd == "TREE-QMC-n0") { mthd2 <- xdf$TQMCn0xFPR }
        else {
                print("ERROR")
                exit()
            }
        }

    # Count ties
    diff <- mthd2 - mthd1  # positive means TQMC-wh-n2 is better
    ntie <- sum(diff == 0)
    nother <- sum(diff < 0)
    ntreeqmc <- sum(diff > 0)

    note <- ""

    if ((ntreeqmc + nother) == 0) {
        pval <- "NA"
        stars <- ""
        mc <- ""
    } else {
        # Run test
        xwsr <-  pvalue(wilcoxsign_test(mthd1 ~ mthd2,
                            zero.method="Wilcoxon",
                            distribution = "exact",
                            paired = TRUE,
                            alternative = "two.sided",
                            conf.int = TRUE))
        ywsr <- pvalue(wilcoxsign_test(mthd1 ~ mthd2, 
                            zero.method="Pratt",
                            distribution = "exact",
                            paired = TRUE,
                            alternative = "two.sided",
                            conf.int = TRUE))

        # Take greater of two tests to be conservative
        if (xwsr < ywsr) {
            wsr <- ywsr
        } else {
            wsr <- xwsr
        }

        stars = ""
        if (wsr < 0.000005) {
            stars = "*****"
        } else if (wsr < 0.00005) {
            stars = "****"
        } else if (wsr < 0.0005) {
            stars = "***"
        } else if (wsr < 0.005) {
            stars = "**"
        } else if (wsr < threshold0) {
            stars = "*"
        }

        mc <- ""
        if (wsr < threshold_bonferoni) {
            mc <- "MC"
        }

        if ((mean(diff) < 0) & (wsr < threshold0)) {
            note <- paste(mthd, " better", sep="")
        }

        pval <- format(wsr, scientific = TRUE, digits=1)
    } # CONDITIONAL

    if ((doit != "varypsiz") &
        (ntax == 50) &
        (ngen == 1000) &
        (nbps == 100) &
        (blsc == 1.0) &
        (psiz == "50000000") &
        (miss == 0.6)) {
        tmp <- note
        note <- paste(tmp, "& (dup)")
        ntest <- 0
    } else {
        ntest <- 1
    }

    writeme <- paste("&", name,
                     "&",
                     "&", as.character(ntreeqmc),
                     "&", as.character(nother),
                     "&", as.character(ntie),
                     "&",
                     "&", pval ,
                     "&", stars,
                     "&", mc,
                     "&", note, "\\\\")

    message(writeme)

    return(ntest)
}


# MAIN FUNCTION
mthds <- c("ASTEROID", "ASTER", "ASTRID",
           "TREE-QMC-n2-shared", "TREE-QMC-n1", "TREE-QMC-n0")
varys <- c("varypsiz", "varyntax", "varyngen",
           "varynbps", "varyblsc", "varymiss")
emets <- c("FNR", "FPR")

ntests <- 312 # CHECK FROM RUNNING CODE
threshold0 <- 0.05
threshold_bonferoni <- threshold0 / ntests

count_tests <- 0

for (mthd in mthds) {
    message("\\begin{table}[!h]")
    message(paste("\\caption[Statistical testing for TREE-QMC-n2 vs", mthd)) 
    message("on Asteroid data]{\\textbf{Testing for differences between TREE-QMC-n2 vs")
    message(paste(mthd, "on the Asteroid data.} BET is the number of replicates for which"))
    message(paste("TREE-QMC-n2 has lower species tree error and thus is better than ", mthd, ", and WOR is the opposite." , sep="")) 
    message("Significance is evaluated using Wilcoxon signed-rank tests on the error rates.")
    message("The symbols *, **, ***, ****, ***** indicate significance at $p <$ 0.5, 0.005, and so on.")
    message(paste("MC indicates significance after Bonferroni correction, i.e., $p < $",
                  as.character(threshold0), "/",
                  as.character(ntests), "=", 
                  format(threshold_bonferoni, scientific = TRUE, digits=1)))
    message(paste("for the ",
                  as.character(ntests),
                  "tests made on the Asteroid data."))
    message("}")
    message("\\centering")
    message("\\scriptsize")
    message("\\begin{tabular}{l r l r r r l c c c l l c}")
    message("\\toprule") 
    message(paste("\\multicolumn{11}{c}{\\textbf{TREE-QMC-n2  vs ", mthd, "}} \\\\"))
    message("\\midrule")
    message("& & & BET & WOR & TIE & & p-val & sig & MC & note \\\\")

    for (emet in emets) {
        message("\\midrule")

        if (emet == "FNR") {
            message("\\multicolumn{11}{c}{\\textbf{False negative error rate}} \\\\")
        } else if (emet == "FPR") {
            message("\\multicolumn{11}{c}{\\textbf{False positive error rate}} \\\\")
        } else {
            print("Unrecognized error metric")
            exit()
        }

        for (vary in varys) {
            message("\\midrule")

            data <- read.csv(paste("../csvs/data-", vary, "-for-testing.csv", sep=""))

            data$NTAX <- as.factor(data$NTAX)
            data$NGEN <- as.factor(data$NGEN)
            data$NBPS <- as.factor(data$NBPS)
            data$BLSC <- as.factor(data$BLSC)
            data$PSIZ <- as.factor(data$PSIZ)
            data$MISS <- as.factor(data$MISS)
            data$REPL <- as.factor(data$REPL)

            ntax <- 50
            ngen <- 1000
            nbps <- 100
            blsc <- 1.0
            psiz <- "50000000"
            miss <- 0.6  # ms always same as mf

            if (vary == "varypsiz") {
                message(paste("\\multirow{ 5}{2cm}{Varying population size}"))
                psizs <- c("10", "50000000", "100000000", "500000000", "1000000000")
                for (psiz in psizs) {
                    #print(psiz)
                    didtest <- run_tests_helper(data, mthd, ntax, ngen, nbps, blsc, psiz, miss,
                                                psiz, emet, vary, threshold_bonferoni)
                    count_tests <- count_tests + didtest
                    #print("DONE")
                }
            } else if (vary == "varyntax") {
                message(paste("\\multirow{ 6}{2cm}{Varying number of taxa}"))
                ntaxs <- c(25, 75, 50, 100, 125, 150)
                for (ntax in ntaxs) {
                    didtest <- run_tests_helper(data, mthd, ntax, ngen, nbps, blsc, psiz, miss,
                                                ntax, emet, vary, threshold_bonferoni)
                    count_tests <- count_tests + didtest
                }
            } else if (vary == "varyngen") {
                message(paste("\\multirow{ 4}{2cm}{Varying number of genes}"))
                ngens <- c(250, 500, 1000, 2000)
                for (ngen in ngens) {
                    didtest <- run_tests_helper(data, mthd, ntax, ngen, nbps, blsc, psiz, miss,
                                                ngen, emet, vary, threshold_bonferoni)
                    count_tests <- count_tests + didtest
                }
            } else if (vary == "varynbps") {
                message(paste("\\multirow{ 4}{2cm}{Varying sequence length}"))
                nbpss <- c(50, 100, 200, 500)
                for (nbps in nbpss) {
                    didtest <- run_tests_helper(data, mthd, ntax, ngen, nbps, blsc, psiz, miss,
                                                nbps, emet, vary, threshold_bonferoni)
                    count_tests <- count_tests + didtest
                }
            } else if (vary == "varyblsc") {
                message(paste("\\multirow{ 6}{2cm}{Varying branch length scaler}"))
                blscs <- c(0.05, 0.1, 1, 10, 100, 200)
                for (blsc in blscs) {
                    didtest <- run_tests_helper(data, mthd, ntax, ngen, nbps, blsc, psiz, miss,
                                                blsc, emet, vary, threshold_bonferoni)
                    count_tests <- count_tests + didtest
                }
            } else if (vary == "varymiss") {
                message(paste("\\multirow{ 6}{2cm}{Varying missingness parameter}"))
                misss <- c(0.5, 0.55, 0.6, 0.65, 0.7, 0.75)
                for (miss in misss) {
                    didtest <- run_tests_helper(data, mthd, ntax, ngen, nbps, blsc, psiz, miss,
                                                miss, emet, vary, threshold_bonferoni)
                    count_tests <- count_tests + didtest
                }
            } else {
                print(paste("Unrecognized experiment", vary))
                exit()
            }
        }  # Varying experiment
    }  # varying emet
    message("\\bottomrule")
    message("\\end{tabular}")
    message("\\end{table}")
    message("")
}  # Varying method

print(paste("did", count_tests, "tests"))
exit()


\begin{table}[!h]
\caption[Statistical testing for TREE-QMC-n2 vs ASTEROID
on Asteroid data]{\textbf{Testing for differences between TREE-QMC-n2 vs
ASTEROID on the Asteroid data.} BET is the number of replicates for which
TREE-QMC-n2 has lower species tree error and thus is better than ASTEROID, and WOR is the opposite.
Significance is evaluated using Wilcoxon signed-rank tests on the error rates.
The symbols *, **, ***, ****, ***** indicate significance at $p <$ 0.5, 0.005, and so on.
MC indicates significance after Bonferroni correction, i.e., $p < $ 0.05 / 312 = 2e-04
for the  312 tests made on the Asteroid data.
}
\centering
\scriptsize
\begin{tabular}{l r l r r r l c c c l l c}
\toprule
\multicolumn{11}{c}{\textbf{TREE-QMC-n2  vs  ASTEROID }} \\
\midrule
& & & BET & WOR & TIE & & p-val & sig & MC & note \\
\midrule
\multicolumn{11}{c}{\textbf{False negative error rate}} \\
\midrule
\multirow{ 5}{2cm}{Varying population size}
& 10 & & 15 & 24 & 11 & & 0.04 & * &  & ASTEROID better \\
& 50000000 & & 12 & 27 & 11 & & 0.01 & * &  & ASTEROID better \\
& 100000000 & & 20 & 20 & 10 & & 0.4 &  &  &  \\
& 500000000 & & 26 & 20 & 4 & & 0.3 &  &  &  \\
& 1000000000 & & 16 & 25 & 9 & & 0.4 &  &  &  \\
\midrule
\multirow{ 6}{2cm}{Varying number of taxa}
& 25 & & 13 & 18 & 19 & & 0.3 &  &  &  \\
& 75 & & 14 & 27 & 9 & & 0.008 & * &  & ASTEROID better \\
& 50 & & 12 & 27 & 11 & & 0.01 & * &  & ASTEROID better & (dup) \\
& 100 & & 13 & 28 & 9 & & 0.1 &  &  &  \\
& 125 & & 13 & 29 & 8 & & 0.003 & ** &  & ASTEROID better \\
& 150 & & 15 & 32 & 3 & & 0.002 & ** &  & ASTEROID better \\
\midrule
\multirow{ 4}{2cm}{Varying number of genes}
& 250 & & 12 & 36 & 2 & & 3e-06 & ***** & MC & ASTEROID better \\
& 500 & & 15 & 28 & 7 & & 0.004 & ** &  & ASTEROID better \\
& 1000 & & 12 & 27 & 11 & & 0.01 & * &  & ASTEROID better & (dup) \\
& 2000 & & 13 & 15 & 22 & & 0.8 &  &  &  \\
\midrule
\multirow{ 4}{2cm}{Varying sequence length}
& 50 & & 13 & 28 & 9 & & 0.02 & * &  & ASTEROID better \\
& 100 & & 12 & 27 & 11 & & 0.01 & * &  & ASTEROID better & (dup) \\
& 200 & & 13 & 24 & 13 & & 0.03 & * &  & ASTEROID better \\
& 500 & & 15 & 20 & 15 & & 0.3 &  &  &  \\
\midrule
\multirow{ 6}{2cm}{Varying branch length scaler}
& 0.05 & & 14 & 26 & 10 & & 0.07 &  &  &  \\
& 0.1 & & 9 & 34 & 7 & & 7e-04 & ** &  & ASTEROID better \\
& 1 & & 12 & 27 & 11 & & 0.01 & * &  & ASTEROID better & (dup) \\
& 10 & & 18 & 18 & 14 & & 0.7 &  &  &  \\
& 100 & & 21 & 25 & 4 & & 0.3 &  &  &  \\
& 200 & & 12 & 29 & 9 & & 0.002 & ** &  & ASTEROID better \\
\midrule
\multirow{ 6}{2cm}{Varying missingness parameter}
& 0.5 & & 14 & 23 & 13 & & 0.3 &  &  &  \\
& 0.55 & & 15 & 17 & 18 & & 0.8 &  &  &  \\
& 0.6 & & 12 & 27 & 11 & & 0.01 & * &  & ASTEROID better & (dup) \\
& 0.65 & & 15 & 24 & 11 & & 0.08 &  &  &  \\
& 0.7 & & 10 & 32 & 8 & & 5e-06 & ***** & MC & ASTEROID better \\
& 0.75 & & 7 & 40 & 3 & & 5e-09 & ***** & MC & ASTEROID better \\
\midrule
\multicolumn{11}{c}{\textbf{False positive error rate}} \\
\midrule
\multirow{ 5}{2cm}{Varying population size}
& 10 & & 20 & 23 & 7 & & 1 &  &  &  \\
& 50000000 & & 24 & 23 & 3 & & 0.7 &  &  &  \\
& 100000000 & & 30 & 19 & 1 & & 0.2 &  &  &  \\
& 500000000 & & 30 & 18 & 2 & & 0.04 & * &  &  \\
& 1000000000 & & 22 & 24 & 4 & & 0.8 &  &  &  \\
\midrule
\multirow{ 6}{2cm}{Varying number of taxa}
& 25 & & 20 & 16 & 14 & & 0.7 &  &  &  \\
& 75 & & 23 & 25 & 2 & & 0.8 &  &  &  \\
& 50 & & 24 & 23 & 3 & & 0.7 &  &  &  & (dup) \\
& 100 & & 29 & 20 & 1 & & 0.03 & * &  &  \\
& 125 & & 30 & 19 & 1 & & 0.01 & * &  &  \\
& 150 & & 31 & 19 & 0 & & 0.1 &  &  &  \\
\midrule
\multirow{ 4}{2cm}{Varying number of genes}
& 250 & & 28 & 22 & 0 & & 0.9 &  &  &  \\
& 500 & & 26 & 23 & 1 & & 0.3 &  &  &  \\
& 1000 & & 24 & 23 & 3 & & 0.7 &  &  &  & (dup) \\
& 2000 & & 21 & 13 & 16 & & 0.09 &  &  &  \\
\midrule
\multirow{ 4}{2cm}{Varying sequence length}
& 50 & & 19 & 27 & 4 & & 0.3 &  &  &  \\
& 100 & & 24 & 23 & 3 & & 0.7 &  &  &  & (dup) \\
& 200 & & 27 & 19 & 4 & & 0.2 &  &  &  \\
& 500 & & 32 & 11 & 7 & & 8e-04 & ** &  &  \\
\midrule
\multirow{ 6}{2cm}{Varying branch length scaler}
& 0.05 & & 26 & 23 & 1 & & 0.9 &  &  &  \\
& 0.1 & & 16 & 29 & 5 & & 0.04 & * &  & ASTEROID better \\
& 1 & & 24 & 23 & 3 & & 0.7 &  &  &  & (dup) \\
& 10 & & 31 & 16 & 3 & & 0.02 & * &  &  \\
& 100 & & 25 & 23 & 2 & & 0.7 &  &  &  \\
& 200 & & 17 & 29 & 4 & & 0.02 & * &  & ASTEROID better \\
\midrule
\multirow{ 6}{2cm}{Varying missingness parameter}
& 0.5 & & 17 & 23 & 10 & & 0.3 &  &  &  \\
& 0.55 & & 22 & 16 & 12 & & 0.3 &  &  &  \\
& 0.6 & & 24 & 23 & 3 & & 0.7 &  &  &  & (dup) \\
& 0.65 & & 39 & 11 & 0 & & 1e-05 & **** & MC &  \\
& 0.7 & & 33 & 17 & 0 & & 0.07 &  &  &  \\
& 0.75 & & 28 & 22 & 0 & & 0.2 &  &  &  \\
\bottomrule
\end{tabular}
\end{table}

\begin{table}[!h]
\caption[Statistical testing for TREE-QMC-n2 vs ASTER
on Asteroid data]{\textbf{Testing for differences between TREE-QMC-n2 vs
ASTER on the Asteroid data.} BET is the number of replicates for which
TREE-QMC-n2 has lower species tree error and thus is better than ASTER, and WOR is the opposite.
Significance is evaluated using Wilcoxon signed-rank tests on the error rates.
The symbols *, **, ***, ****, ***** indicate significance at $p <$ 0.5, 0.005, and so on.
MC indicates significance after Bonferroni correction, i.e., $p < $ 0.05 / 312 = 2e-04
for the  312 tests made on the Asteroid data.
}
\centering
\scriptsize
\begin{tabular}{l r l r r r l c c c l l c}
\toprule
\multicolumn{11}{c}{\textbf{TREE-QMC-n2  vs  ASTER }} \\
\midrule
& & & BET & WOR & TIE & & p-val & sig & MC & note \\
\midrule
\multicolumn{11}{c}{\textbf{False negative error rate}} \\
\midrule
\multirow{ 5}{2cm}{Varying population size}
& 10 & & 30 & 13 & 7 & & 0.005 & * &  &  \\
& 50000000 & & 24 & 19 & 7 & & 0.3 &  &  &  \\
& 100000000 & & 26 & 16 & 8 & & 0.02 & * &  &  \\
& 500000000 & & 32 & 12 & 6 & & 5e-04 & *** &  &  \\
& 1000000000 & & 32 & 16 & 2 & & 0.005 & * &  &  \\
\midrule
\multirow{ 6}{2cm}{Varying number of taxa}
& 25 & & 29 & 9 & 12 & & 8e-04 & ** &  &  \\
& 75 & & 26 & 18 & 6 & & 0.1 &  &  &  \\
& 50 & & 24 & 19 & 7 & & 0.3 &  &  &  & (dup) \\
& 100 & & 30 & 15 & 5 & & 0.01 & * &  &  \\
& 125 & & 31 & 16 & 3 & & 0.003 & ** &  &  \\
& 150 & & 38 & 6 & 6 & & 3e-06 & ***** & MC &  \\
\midrule
\multirow{ 4}{2cm}{Varying number of genes}
& 250 & & 28 & 16 & 6 & & 0.1 &  &  &  \\
& 500 & & 30 & 14 & 6 & & 0.006 & * &  &  \\
& 1000 & & 24 & 19 & 7 & & 0.3 &  &  &  & (dup) \\
& 2000 & & 34 & 9 & 7 & & 0.001 & ** &  &  \\
\midrule
\multirow{ 4}{2cm}{Varying sequence length}
& 50 & & 29 & 11 & 10 & & 0.001 & ** &  &  \\
& 100 & & 24 & 19 & 7 & & 0.3 &  &  &  & (dup) \\
& 200 & & 17 & 26 & 7 & & 0.4 &  &  &  \\
& 500 & & 21 & 17 & 12 & & 0.3 &  &  &  \\
\midrule
\multirow{ 6}{2cm}{Varying branch length scaler}
& 0.05 & & 27 & 17 & 6 & & 0.09 &  &  &  \\
& 0.1 & & 27 & 15 & 8 & & 0.05 &  &  &  \\
& 1 & & 24 & 19 & 7 & & 0.3 &  &  &  & (dup) \\
& 10 & & 29 & 13 & 8 & & 6e-04 & ** &  &  \\
& 100 & & 37 & 10 & 3 & & 2e-06 & ***** & MC &  \\
& 200 & & 25 & 21 & 4 & & 0.6 &  &  &  \\
\midrule
\multirow{ 6}{2cm}{Varying missingness parameter}
& 0.5 & & 27 & 11 & 12 & & 0.02 & * &  &  \\
& 0.55 & & 31 & 9 & 10 & & 5e-04 & *** &  &  \\
& 0.6 & & 24 & 19 & 7 & & 0.3 &  &  &  & (dup) \\
& 0.65 & & 29 & 7 & 14 & & 2e-04 & *** & MC &  \\
& 0.7 & & 22 & 20 & 8 & & 0.4 &  &  &  \\
& 0.75 & & 28 & 15 & 7 & & 0.01 & * &  &  \\
\midrule
\multicolumn{11}{c}{\textbf{False positive error rate}} \\
\midrule
\multirow{ 5}{2cm}{Varying population size}
& 10 & & 36 & 12 & 2 & & 3e-05 & **** & MC &  \\
& 50000000 & & 28 & 18 & 4 & & 0.009 & * &  &  \\
& 100000000 & & 36 & 14 & 0 & & 5e-05 & *** & MC &  \\
& 500000000 & & 40 & 10 & 0 & & 1e-05 & **** & MC &  \\
& 1000000000 & & 32 & 16 & 2 & & 2e-04 & *** &  &  \\
\midrule
\multirow{ 6}{2cm}{Varying number of taxa}
& 25 & & 36 & 7 & 7 & & 8e-06 & **** & MC &  \\
& 75 & & 33 & 16 & 1 & & 3e-04 & *** &  &  \\
& 50 & & 28 & 18 & 4 & & 0.009 & * &  &  & (dup) \\
& 100 & & 39 & 10 & 1 & & 4e-08 & ***** & MC &  \\
& 125 & & 41 & 9 & 0 & & 2e-08 & ***** & MC &  \\
& 150 & & 46 & 4 & 0 & & 4e-13 & ***** & MC &  \\
\midrule
\multirow{ 4}{2cm}{Varying number of genes}
& 250 & & 40 & 10 & 0 & & 2e-08 & ***** & MC &  \\
& 500 & & 39 & 11 & 0 & & 3e-07 & ***** & MC &  \\
& 1000 & & 28 & 18 & 4 & & 0.009 & * &  &  & (dup) \\
& 2000 & & 39 & 8 & 3 & & 2e-06 & ***** & MC &  \\
\midrule
\multirow{ 4}{2cm}{Varying sequence length}
& 50 & & 36 & 10 & 4 & & 2e-05 & **** & MC &  \\
& 100 & & 28 & 18 & 4 & & 0.009 & * &  &  & (dup) \\
& 200 & & 29 & 18 & 3 & & 0.009 & * &  &  \\
& 500 & & 33 & 14 & 3 & & 2e-04 & *** &  &  \\
\midrule
\multirow{ 6}{2cm}{Varying branch length scaler}
& 0.05 & & 32 & 16 & 2 & & 0.002 & ** &  &  \\
& 0.1 & & 29 & 15 & 6 & & 0.002 & ** &  &  \\
& 1 & & 28 & 18 & 4 & & 0.009 & * &  &  & (dup) \\
& 10 & & 37 & 12 & 1 & & 5e-07 & ***** & MC &  \\
& 100 & & 41 & 9 & 0 & & 8e-08 & ***** & MC &  \\
& 200 & & 28 & 20 & 2 & & 0.2 &  &  &  \\
\midrule
\multirow{ 6}{2cm}{Varying missingness parameter}
& 0.5 & & 28 & 11 & 11 & & 0.009 & * &  &  \\
& 0.55 & & 36 & 9 & 5 & & 3e-05 & **** & MC &  \\
& 0.6 & & 28 & 18 & 4 & & 0.009 & * &  &  & (dup) \\
& 0.65 & & 45 & 4 & 1 & & 1e-12 & ***** & MC &  \\
& 0.7 & & 43 & 7 & 0 & & 8e-09 & ***** & MC &  \\
& 0.75 & & 45 & 5 & 0 & & 3e-10 & ***** & MC &  \\
\bottomrule
\end{tabular}
\end{table}

\begin{table}[!h]
\caption[Statistical testing for TREE-QMC-n2 vs ASTRID
on Asteroid data]{\textbf{Testing for differences between TREE-QMC-n2 vs
ASTRID on the Asteroid data.} BET is the number of replicates for which
TREE-QMC-n2 has lower species tree error and thus is better than ASTRID, and WOR is the opposite.
Significance is evaluated using Wilcoxon signed-rank tests on the error rates.
The symbols *, **, ***, ****, ***** indicate significance at $p <$ 0.5, 0.005, and so on.
MC indicates significance after Bonferroni correction, i.e., $p < $ 0.05 / 312 = 2e-04
for the  312 tests made on the Asteroid data.
}
\centering
\scriptsize
\begin{tabular}{l r l r r r l c c c l l c}
\toprule
\multicolumn{11}{c}{\textbf{TREE-QMC-n2  vs  ASTRID }} \\
\midrule
& & & BET & WOR & TIE & & p-val & sig & MC & note \\
\midrule
\multicolumn{11}{c}{\textbf{False negative error rate}} \\
\midrule
\multirow{ 5}{2cm}{Varying population size}
& 10 & & 40 & 5 & 5 & & 3e-09 & ***** & MC &  \\
& 50000000 & & 35 & 10 & 5 & & 1e-06 & ***** & MC &  \\
& 100000000 & & 38 & 8 & 4 & & 1e-06 & ***** & MC &  \\
& 500000000 & & 34 & 11 & 5 & & 5e-05 & **** & MC &  \\
& 1000000000 & & 30 & 15 & 5 & & 0.05 &  &  &  \\
\midrule
\multirow{ 6}{2cm}{Varying number of taxa}
& 25 & & 25 & 14 & 11 & & 0.04 & * &  &  \\
& 75 & & 36 & 9 & 5 & & 5e-07 & ***** & MC &  \\
& 50 & & 35 & 10 & 5 & & 1e-06 & ***** & MC &  & (dup) \\
& 100 & & 40 & 8 & 2 & & 6e-09 & ***** & MC &  \\
& 125 & & 45 & 4 & 1 & & 3e-11 & ***** & MC &  \\
& 150 & & 48 & 2 & 0 & & 5e-13 & ***** & MC &  \\
\midrule
\multirow{ 4}{2cm}{Varying number of genes}
& 250 & & 38 & 9 & 3 & & 1e-05 & **** & MC &  \\
& 500 & & 40 & 7 & 3 & & 6e-08 & ***** & MC &  \\
& 1000 & & 35 & 10 & 5 & & 1e-06 & ***** & MC &  & (dup) \\
& 2000 & & 36 & 8 & 6 & & 3e-07 & ***** & MC &  \\
\midrule
\multirow{ 4}{2cm}{Varying sequence length}
& 50 & & 30 & 17 & 3 & & 0.007 & * &  &  \\
& 100 & & 35 & 10 & 5 & & 1e-06 & ***** & MC &  & (dup) \\
& 200 & & 41 & 6 & 3 & & 9e-08 & ***** & MC &  \\
& 500 & & 42 & 4 & 4 & & 6e-10 & ***** & MC &  \\
\midrule
\multirow{ 6}{2cm}{Varying branch length scaler}
& 0.05 & & 30 & 14 & 6 & & 0.005 & * &  &  \\
& 0.1 & & 35 & 10 & 5 & & 1e-05 & **** & MC &  \\
& 1 & & 35 & 10 & 5 & & 1e-06 & ***** & MC &  & (dup) \\
& 10 & & 34 & 13 & 3 & & 3e-05 & **** & MC &  \\
& 100 & & 34 & 14 & 2 & & 7e-04 & ** &  &  \\
& 200 & & 26 & 23 & 1 & & 0.6 &  &  &  \\
\midrule
\multirow{ 6}{2cm}{Varying missingness parameter}
& 0.5 & & 29 & 10 & 11 & & 4e-04 & *** &  &  \\
& 0.55 & & 32 & 15 & 3 & & 0.003 & ** &  &  \\
& 0.6 & & 35 & 10 & 5 & & 1e-06 & ***** & MC &  & (dup) \\
& 0.65 & & 35 & 9 & 6 & & 4e-08 & ***** & MC &  \\
& 0.7 & & 43 & 3 & 4 & & 2e-11 & ***** & MC &  \\
& 0.75 & & 46 & 1 & 3 & & 7e-14 & ***** & MC &  \\
\midrule
\multicolumn{11}{c}{\textbf{False positive error rate}} \\
\midrule
\multirow{ 5}{2cm}{Varying population size}
& 10 & & 44 & 4 & 2 & & 7e-12 & ***** & MC &  \\
& 50000000 & & 39 & 9 & 2 & & 1e-08 & ***** & MC &  \\
& 100000000 & & 40 & 8 & 2 & & 3e-08 & ***** & MC &  \\
& 500000000 & & 38 & 10 & 2 & & 5e-06 & **** & MC &  \\
& 1000000000 & & 34 & 15 & 1 & & 0.006 & * &  &  \\
\midrule
\multirow{ 6}{2cm}{Varying number of taxa}
& 25 & & 31 & 12 & 7 & & 0.002 & ** &  &  \\
& 75 & & 42 & 7 & 1 & & 9e-10 & ***** & MC &  \\
& 50 & & 39 & 9 & 2 & & 1e-08 & ***** & MC &  & (dup) \\
& 100 & & 43 & 6 & 1 & & 1e-11 & ***** & MC &  \\
& 125 & & 47 & 3 & 0 & & 2e-13 & ***** & MC &  \\
& 150 & & 49 & 1 & 0 & & 9e-15 & ***** & MC &  \\
\midrule
\multirow{ 4}{2cm}{Varying number of genes}
& 250 & & 45 & 5 & 0 & & 1e-10 & ***** & MC &  \\
& 500 & & 43 & 7 & 0 & & 1e-10 & ***** & MC &  \\
& 1000 & & 39 & 9 & 2 & & 1e-08 & ***** & MC &  & (dup) \\
& 2000 & & 40 & 6 & 4 & & 1e-09 & ***** & MC &  \\
\midrule
\multirow{ 4}{2cm}{Varying sequence length}
& 50 & & 31 & 17 & 2 & & 3e-04 & *** &  &  \\
& 100 & & 39 & 9 & 2 & & 1e-08 & ***** & MC &  & (dup) \\
& 200 & & 45 & 5 & 0 & & 1e-10 & ***** & MC &  \\
& 500 & & 45 & 2 & 3 & & 4e-12 & ***** & MC &  \\
\midrule
\multirow{ 6}{2cm}{Varying branch length scaler}
& 0.05 & & 35 & 14 & 1 & & 3e-04 & *** &  &  \\
& 0.1 & & 38 & 10 & 2 & & 2e-07 & ***** & MC &  \\
& 1 & & 39 & 9 & 2 & & 1e-08 & ***** & MC &  & (dup) \\
& 10 & & 39 & 10 & 1 & & 1e-07 & ***** & MC &  \\
& 100 & & 35 & 14 & 1 & & 2e-04 & *** &  &  \\
& 200 & & 27 & 23 & 0 & & 0.3 &  &  &  \\
\midrule
\multirow{ 6}{2cm}{Varying missingness parameter}
& 0.5 & & 31 & 10 & 9 & & 2e-04 & *** &  &  \\
& 0.55 & & 33 & 15 & 2 & & 9e-04 & ** &  &  \\
& 0.6 & & 39 & 9 & 2 & & 1e-08 & ***** & MC &  & (dup) \\
& 0.65 & & 45 & 5 & 0 & & 4e-13 & ***** & MC &  \\
& 0.7 & & 49 & 1 & 0 & & 4e-15 & ***** & MC &  \\
& 0.75 & & 49 & 1 & 0 & & 4e-15 & ***** & MC &  \\
\bottomrule
\end{tabular}
\end{table}

\begin{table}[!h]
\caption[Statistical testing for TREE-QMC-n2 vs TREE-QMC-n2-shared
on Asteroid data]{\textbf{Testing for differences between TREE-QMC-n2 vs
TREE-QMC-n2-shared on the Asteroid data.} BET is the number of replicates for which
TREE-QMC-n2 has lower species tree error and thus is better than TREE-QMC-n2-shared, and WOR is the opposite.
Significance is evaluated using Wilcoxon signed-rank tests on the error rates.
The symbols *, **, ***, ****, ***** indicate significance at $p <$ 0.5, 0.005, and so on.
MC indicates significance after Bonferroni correction, i.e., $p < $ 0.05 / 312 = 2e-04
for the  312 tests made on the Asteroid data.
}
\centering
\scriptsize
\begin{tabular}{l r l r r r l c c c l l c}
\toprule
\multicolumn{11}{c}{\textbf{TREE-QMC-n2  vs  TREE-QMC-n2-shared }} \\
\midrule
& & & BET & WOR & TIE & & p-val & sig & MC & note \\
\midrule
\multicolumn{11}{c}{\textbf{False negative error rate}} \\
\midrule
\multirow{ 5}{2cm}{Varying population size}
& 10 & & 41 & 5 & 4 & & 1e-09 & ***** & MC &  \\
& 50000000 & & 37 & 7 & 6 & & 5e-08 & ***** & MC &  \\
& 100000000 & & 45 & 1 & 4 & & 2e-13 & ***** & MC &  \\
& 500000000 & & 44 & 3 & 3 & & 2e-10 & ***** & MC &  \\
& 1000000000 & & 46 & 4 & 0 & & 4e-11 & ***** & MC &  \\
\midrule
\multirow{ 6}{2cm}{Varying number of taxa}
& 25 & & 28 & 9 & 13 & & 1e-04 & *** & MC &  \\
& 75 & & 46 & 2 & 2 & & 1e-12 & ***** & MC &  \\
& 50 & & 37 & 7 & 6 & & 5e-08 & ***** & MC &  & (dup) \\
& 100 & & 50 & 0 & 0 & & 2e-15 & ***** & MC &  \\
& 125 & & 50 & 0 & 0 & & 2e-15 & ***** & MC &  \\
& 150 & & 50 & 0 & 0 & & 2e-15 & ***** & MC &  \\
\midrule
\multirow{ 4}{2cm}{Varying number of genes}
& 250 & & 41 & 5 & 4 & & 4e-10 & ***** & MC &  \\
& 500 & & 41 & 6 & 3 & & 3e-10 & ***** & MC &  \\
& 1000 & & 37 & 7 & 6 & & 5e-08 & ***** & MC &  & (dup) \\
& 2000 & & 41 & 4 & 5 & & 8e-11 & ***** & MC &  \\
\midrule
\multirow{ 4}{2cm}{Varying sequence length}
& 50 & & 41 & 3 & 6 & & 2e-11 & ***** & MC &  \\
& 100 & & 37 & 7 & 6 & & 5e-08 & ***** & MC &  & (dup) \\
& 200 & & 37 & 7 & 6 & & 9e-08 & ***** & MC &  \\
& 500 & & 40 & 4 & 6 & & 4e-08 & ***** & MC &  \\
\midrule
\multirow{ 6}{2cm}{Varying branch length scaler}
& 0.05 & & 49 & 1 & 0 & & 5e-15 & ***** & MC &  \\
& 0.1 & & 48 & 2 & 0 & & 1e-13 & ***** & MC &  \\
& 1 & & 37 & 7 & 6 & & 5e-08 & ***** & MC &  & (dup) \\
& 10 & & 43 & 6 & 1 & & 3e-10 & ***** & MC &  \\
& 100 & & 48 & 2 & 0 & & 2e-14 & ***** & MC &  \\
& 200 & & 46 & 4 & 0 & & 1e-11 & ***** & MC &  \\
\midrule
\multirow{ 6}{2cm}{Varying missingness parameter}
& 0.5 & & 34 & 7 & 9 & & 2e-06 & ***** & MC &  \\
& 0.55 & & 39 & 5 & 6 & & 1e-09 & ***** & MC &  \\
& 0.6 & & 37 & 7 & 6 & & 5e-08 & ***** & MC &  & (dup) \\
& 0.65 & & 45 & 5 & 0 & & 1e-10 & ***** & MC &  \\
& 0.7 & & 44 & 1 & 5 & & 2e-13 & ***** & MC &  \\
& 0.75 & & 41 & 6 & 3 & & 4e-09 & ***** & MC &  \\
\midrule
\multicolumn{11}{c}{\textbf{False positive error rate}} \\
\midrule
\multirow{ 5}{2cm}{Varying population size}
& 10 & & 42 & 5 & 3 & & 7e-10 & ***** & MC &  \\
& 50000000 & & 39 & 7 & 4 & & 1e-08 & ***** & MC &  \\
& 100000000 & & 46 & 1 & 3 & & 4e-14 & ***** & MC &  \\
& 500000000 & & 45 & 4 & 1 & & 2e-10 & ***** & MC &  \\
& 1000000000 & & 46 & 4 & 0 & & 3e-11 & ***** & MC &  \\
\midrule
\multirow{ 6}{2cm}{Varying number of taxa}
& 25 & & 28 & 11 & 11 & & 2e-04 & *** &  &  \\
& 75 & & 46 & 3 & 1 & & 7e-13 & ***** & MC &  \\
& 50 & & 39 & 7 & 4 & & 1e-08 & ***** & MC &  & (dup) \\
& 100 & & 50 & 0 & 0 & & 2e-15 & ***** & MC &  \\
& 125 & & 50 & 0 & 0 & & 2e-15 & ***** & MC &  \\
& 150 & & 50 & 0 & 0 & & 2e-15 & ***** & MC &  \\
\midrule
\multirow{ 4}{2cm}{Varying number of genes}
& 250 & & 41 & 8 & 1 & & 8e-10 & ***** & MC &  \\
& 500 & & 45 & 5 & 0 & & 9e-11 & ***** & MC &  \\
& 1000 & & 39 & 7 & 4 & & 1e-08 & ***** & MC &  & (dup) \\
& 2000 & & 42 & 4 & 4 & & 3e-10 & ***** & MC &  \\
\midrule
\multirow{ 4}{2cm}{Varying sequence length}
& 50 & & 42 & 6 & 2 & & 1e-11 & ***** & MC &  \\
& 100 & & 39 & 7 & 4 & & 1e-08 & ***** & MC &  & (dup) \\
& 200 & & 38 & 7 & 5 & & 1e-07 & ***** & MC &  \\
& 500 & & 40 & 6 & 4 & & 1e-08 & ***** & MC &  \\
\midrule
\multirow{ 6}{2cm}{Varying branch length scaler}
& 0.05 & & 49 & 1 & 0 & & 4e-15 & ***** & MC &  \\
& 0.1 & & 48 & 2 & 0 & & 6e-14 & ***** & MC &  \\
& 1 & & 39 & 7 & 4 & & 1e-08 & ***** & MC &  & (dup) \\
& 10 & & 44 & 6 & 0 & & 6e-11 & ***** & MC &  \\
& 100 & & 49 & 1 & 0 & & 1e-14 & ***** & MC &  \\
& 200 & & 47 & 3 & 0 & & 9e-12 & ***** & MC &  \\
\midrule
\multirow{ 6}{2cm}{Varying missingness parameter}
& 0.5 & & 34 & 7 & 9 & & 3e-06 & ***** & MC &  \\
& 0.55 & & 38 & 6 & 6 & & 2e-09 & ***** & MC &  \\
& 0.6 & & 39 & 7 & 4 & & 1e-08 & ***** & MC &  & (dup) \\
& 0.65 & & 46 & 4 & 0 & & 1e-11 & ***** & MC &  \\
& 0.7 & & 46 & 3 & 1 & & 5e-13 & ***** & MC &  \\
& 0.75 & & 44 & 6 & 0 & & 2e-10 & ***** & MC &  \\
\bottomrule
\end{tabular}
\end{table}

\begin{table}[!h]
\caption[Statistical testing for TREE-QMC-n2 vs TREE-QMC-n1
on Asteroid data]{\textbf{Testing for differences between TREE-QMC-n2 vs
TREE-QMC-n1 on the Asteroid data.} BET is the number of replicates for which
TREE-QMC-n2 has lower species tree error and thus is better than TREE-QMC-n1, and WOR is the opposite.
Significance is evaluated using Wilcoxon signed-rank tests on the error rates.
The symbols *, **, ***, ****, ***** indicate significance at $p <$ 0.5, 0.005, and so on.
MC indicates significance after Bonferroni correction, i.e., $p < $ 0.05 / 312 = 2e-04
for the  312 tests made on the Asteroid data.
}
\centering
\scriptsize
\begin{tabular}{l r l r r r l c c c l l c}
\toprule
\multicolumn{11}{c}{\textbf{TREE-QMC-n2  vs  TREE-QMC-n1 }} \\
\midrule
& & & BET & WOR & TIE & & p-val & sig & MC & note \\
\midrule
\multicolumn{11}{c}{\textbf{False negative error rate}} \\
\midrule
\multirow{ 5}{2cm}{Varying population size}
& 10 & & 17 & 12 & 21 & & 0.6 &  &  &  \\
& 50000000 & & 11 & 13 & 26 & & 0.9 &  &  &  \\
& 100000000 & & 18 & 16 & 16 & & 0.6 &  &  &  \\
& 500000000 & & 19 & 17 & 14 & & 0.5 &  &  &  \\
& 1000000000 & & 14 & 15 & 21 & & 0.8 &  &  &  \\
\midrule
\multirow{ 6}{2cm}{Varying number of taxa}
& 25 & & 7 & 8 & 35 & & 1 &  &  &  \\
& 75 & & 22 & 15 & 13 & & 0.2 &  &  &  \\
& 50 & & 11 & 13 & 26 & & 0.9 &  &  &  & (dup) \\
& 100 & & 21 & 17 & 12 & & 0.2 &  &  &  \\
& 125 & & 26 & 16 & 8 & & 0.04 & * &  &  \\
& 150 & & 27 & 16 & 7 & & 0.03 & * &  &  \\
\midrule
\multirow{ 4}{2cm}{Varying number of genes}
& 250 & & 18 & 14 & 18 & & 0.5 &  &  &  \\
& 500 & & 20 & 13 & 17 & & 0.2 &  &  &  \\
& 1000 & & 11 & 13 & 26 & & 0.9 &  &  &  & (dup) \\
& 2000 & & 14 & 8 & 28 & & 0.2 &  &  &  \\
\midrule
\multirow{ 4}{2cm}{Varying sequence length}
& 50 & & 16 & 17 & 17 & & 0.5 &  &  &  \\
& 100 & & 11 & 13 & 26 & & 0.9 &  &  &  & (dup) \\
& 200 & & 14 & 14 & 22 & & 1 &  &  &  \\
& 500 & & 6 & 13 & 31 & & 0.09 &  &  &  \\
\midrule
\multirow{ 6}{2cm}{Varying branch length scaler}
& 0.05 & & 15 & 24 & 11 & & 0.3 &  &  &  \\
& 0.1 & & 14 & 22 & 14 & & 0.2 &  &  &  \\
& 1 & & 11 & 13 & 26 & & 0.9 &  &  &  & (dup) \\
& 10 & & 21 & 10 & 19 & & 0.01 & * &  &  \\
& 100 & & 30 & 12 & 8 & & 0.002 & ** &  &  \\
& 200 & & 22 & 19 & 9 & & 0.6 &  &  &  \\
\midrule
\multirow{ 6}{2cm}{Varying missingness parameter}
& 0.5 & & 16 & 12 & 22 & & 0.4 &  &  &  \\
& 0.55 & & 16 & 12 & 22 & & 0.7 &  &  &  \\
& 0.6 & & 11 & 13 & 26 & & 0.9 &  &  &  & (dup) \\
& 0.65 & & 21 & 12 & 17 & & 0.09 &  &  &  \\
& 0.7 & & 21 & 15 & 14 & & 0.3 &  &  &  \\
& 0.75 & & 14 & 22 & 14 & & 0.2 &  &  &  \\
\midrule
\multicolumn{11}{c}{\textbf{False positive error rate}} \\
\midrule
\multirow{ 5}{2cm}{Varying population size}
& 10 & & 17 & 13 & 20 & & 0.5 &  &  &  \\
& 50000000 & & 13 & 13 & 24 & & 1 &  &  &  \\
& 100000000 & & 18 & 17 & 15 & & 0.7 &  &  &  \\
& 500000000 & & 19 & 17 & 14 & & 0.5 &  &  &  \\
& 1000000000 & & 17 & 18 & 15 & & 0.9 &  &  &  \\
\midrule
\multirow{ 6}{2cm}{Varying number of taxa}
& 25 & & 7 & 8 & 35 & & 1 &  &  &  \\
& 75 & & 22 & 16 & 12 & & 0.3 &  &  &  \\
& 50 & & 13 & 13 & 24 & & 1 &  &  &  & (dup) \\
& 100 & & 21 & 17 & 12 & & 0.2 &  &  &  \\
& 125 & & 27 & 16 & 7 & & 0.2 &  &  &  \\
& 150 & & 28 & 16 & 6 & & 0.03 & * &  &  \\
\midrule
\multirow{ 4}{2cm}{Varying number of genes}
& 250 & & 16 & 19 & 15 & & 1 &  &  &  \\
& 500 & & 22 & 16 & 12 & & 0.2 &  &  &  \\
& 1000 & & 13 & 13 & 24 & & 1 &  &  &  & (dup) \\
& 2000 & & 14 & 8 & 28 & & 0.2 &  &  &  \\
\midrule
\multirow{ 4}{2cm}{Varying sequence length}
& 50 & & 16 & 17 & 17 & & 0.7 &  &  &  \\
& 100 & & 13 & 13 & 24 & & 1 &  &  &  & (dup) \\
& 200 & & 14 & 14 & 22 & & 0.8 &  &  &  \\
& 500 & & 6 & 13 & 31 & & 0.06 &  &  &  \\
\midrule
\multirow{ 6}{2cm}{Varying branch length scaler}
& 0.05 & & 16 & 25 & 9 & & 0.4 &  &  &  \\
& 0.1 & & 14 & 23 & 13 & & 0.2 &  &  &  \\
& 1 & & 13 & 13 & 24 & & 1 &  &  &  & (dup) \\
& 10 & & 24 & 10 & 16 & & 0.02 & * &  &  \\
& 100 & & 31 & 12 & 7 & & 0.01 & * &  &  \\
& 200 & & 22 & 19 & 9 & & 0.5 &  &  &  \\
\midrule
\multirow{ 6}{2cm}{Varying missingness parameter}
& 0.5 & & 16 & 12 & 22 & & 0.5 &  &  &  \\
& 0.55 & & 17 & 12 & 21 & & 0.8 &  &  &  \\
& 0.6 & & 13 & 13 & 24 & & 1 &  &  &  & (dup) \\
& 0.65 & & 22 & 12 & 16 & & 0.05 &  &  &  \\
& 0.7 & & 24 & 15 & 11 & & 0.3 &  &  &  \\
& 0.75 & & 18 & 24 & 8 & & 0.3 &  &  &  \\
\bottomrule
\end{tabular}
\end{table}

\begin{table}[!h]
\caption[Statistical testing for TREE-QMC-n2 vs TREE-QMC-n0
on Asteroid data]{\textbf{Testing for differences between TREE-QMC-n2 vs
TREE-QMC-n0 on the Asteroid data.} BET is the number of replicates for which
TREE-QMC-n2 has lower species tree error and thus is better than TREE-QMC-n0, and WOR is the opposite.
Significance is evaluated using Wilcoxon signed-rank tests on the error rates.
The symbols *, **, ***, ****, ***** indicate significance at $p <$ 0.5, 0.005, and so on.
MC indicates significance after Bonferroni correction, i.e., $p < $ 0.05 / 312 = 2e-04
for the  312 tests made on the Asteroid data.
}
\centering
\scriptsize
\begin{tabular}{l r l r r r l c c c l l c}
\toprule
\multicolumn{11}{c}{\textbf{TREE-QMC-n2  vs  TREE-QMC-n0 }} \\
\midrule
& & & BET & WOR & TIE & & p-val & sig & MC & note \\
\midrule
\multicolumn{11}{c}{\textbf{False negative error rate}} \\
\midrule
\multirow{ 5}{2cm}{Varying population size}
& 10 & & 33 & 10 & 7 & & 4e-05 & **** & MC &  \\
& 50000000 & & 30 & 12 & 8 & & 6e-04 & ** &  &  \\
& 100000000 & & 33 & 11 & 6 & & 3e-05 & **** & MC &  \\
& 500000000 & & 33 & 8 & 9 & & 1e-04 & *** & MC &  \\
& 1000000000 & & 35 & 11 & 4 & & 3e-06 & ***** & MC &  \\
\midrule
\multirow{ 6}{2cm}{Varying number of taxa}
& 25 & & 29 & 8 & 13 & & 6e-05 & *** & MC &  \\
& 75 & & 30 & 11 & 9 & & 3e-04 & *** &  &  \\
& 50 & & 30 & 12 & 8 & & 6e-04 & ** &  &  & (dup) \\
& 100 & & 39 & 7 & 4 & & 7e-09 & ***** & MC &  \\
& 125 & & 39 & 6 & 5 & & 1e-07 & ***** & MC &  \\
& 150 & & 46 & 3 & 1 & & 7e-13 & ***** & MC &  \\
\midrule
\multirow{ 4}{2cm}{Varying number of genes}
& 250 & & 34 & 8 & 8 & & 2e-06 & ***** & MC &  \\
& 500 & & 33 & 10 & 7 & & 1e-05 & **** & MC &  \\
& 1000 & & 30 & 12 & 8 & & 6e-04 & ** &  &  & (dup) \\
& 2000 & & 36 & 8 & 6 & & 2e-05 & **** & MC &  \\
\midrule
\multirow{ 4}{2cm}{Varying sequence length}
& 50 & & 31 & 9 & 10 & & 1e-04 & *** & MC &  \\
& 100 & & 30 & 12 & 8 & & 6e-04 & ** &  &  & (dup) \\
& 200 & & 26 & 11 & 13 & & 0.004 & ** &  &  \\
& 500 & & 27 & 8 & 15 & & 0.002 & ** &  &  \\
\midrule
\multirow{ 6}{2cm}{Varying branch length scaler}
& 0.05 & & 29 & 14 & 7 & & 0.01 & * &  &  \\
& 0.1 & & 34 & 10 & 6 & & 7e-05 & *** & MC &  \\
& 1 & & 30 & 12 & 8 & & 6e-04 & ** &  &  & (dup) \\
& 10 & & 34 & 12 & 4 & & 1e-05 & **** & MC &  \\
& 100 & & 37 & 7 & 6 & & 3e-08 & ***** & MC &  \\
& 200 & & 31 & 13 & 6 & & 6e-04 & ** &  &  \\
\midrule
\multirow{ 6}{2cm}{Varying missingness parameter}
& 0.5 & & 33 & 8 & 9 & & 0.001 & ** &  &  \\
& 0.55 & & 31 & 6 & 13 & & 3e-05 & **** & MC &  \\
& 0.6 & & 30 & 12 & 8 & & 6e-04 & ** &  &  & (dup) \\
& 0.65 & & 38 & 8 & 4 & & 3e-07 & ***** & MC &  \\
& 0.7 & & 31 & 11 & 8 & & 0.003 & ** &  &  \\
& 0.75 & & 30 & 15 & 5 & & 0.01 & * &  &  \\
\midrule
\multicolumn{11}{c}{\textbf{False positive error rate}} \\
\midrule
\multirow{ 5}{2cm}{Varying population size}
& 10 & & 35 & 9 & 6 & & 3e-05 & **** & MC &  \\
& 50000000 & & 32 & 13 & 5 & & 5e-04 & ** &  &  \\
& 100000000 & & 35 & 11 & 4 & & 1e-05 & **** & MC &  \\
& 500000000 & & 34 & 11 & 5 & & 1e-04 & *** & MC &  \\
& 1000000000 & & 36 & 14 & 0 & & 2e-06 & ***** & MC &  \\
\midrule
\multirow{ 6}{2cm}{Varying number of taxa}
& 25 & & 30 & 10 & 10 & & 3e-04 & *** &  &  \\
& 75 & & 30 & 13 & 7 & & 4e-04 & *** &  &  \\
& 50 & & 32 & 13 & 5 & & 5e-04 & ** &  &  & (dup) \\
& 100 & & 40 & 8 & 2 & & 3e-09 & ***** & MC &  \\
& 125 & & 42 & 5 & 3 & & 5e-08 & ***** & MC &  \\
& 150 & & 47 & 3 & 0 & & 8e-14 & ***** & MC &  \\
\midrule
\multirow{ 4}{2cm}{Varying number of genes}
& 250 & & 35 & 11 & 4 & & 3e-06 & ***** & MC &  \\
& 500 & & 36 & 9 & 5 & & 4e-06 & ***** & MC &  \\
& 1000 & & 32 & 13 & 5 & & 5e-04 & ** &  &  & (dup) \\
& 2000 & & 36 & 8 & 6 & & 1e-05 & **** & MC &  \\
\midrule
\multirow{ 4}{2cm}{Varying sequence length}
& 50 & & 32 & 9 & 9 & & 5e-05 & **** & MC &  \\
& 100 & & 32 & 13 & 5 & & 5e-04 & ** &  &  & (dup) \\
& 200 & & 28 & 13 & 9 & & 0.004 & ** &  &  \\
& 500 & & 28 & 9 & 13 & & 4e-04 & *** &  &  \\
\midrule
\multirow{ 6}{2cm}{Varying branch length scaler}
& 0.05 & & 33 & 15 & 2 & & 0.004 & ** &  &  \\
& 0.1 & & 35 & 12 & 3 & & 3e-05 & **** & MC &  \\
& 1 & & 32 & 13 & 5 & & 5e-04 & ** &  &  & (dup) \\
& 10 & & 37 & 12 & 1 & & 5e-06 & ***** & MC &  \\
& 100 & & 40 & 8 & 2 & & 3e-08 & ***** & MC &  \\
& 200 & & 35 & 13 & 2 & & 3e-04 & *** &  &  \\
\midrule
\multirow{ 6}{2cm}{Varying missingness parameter}
& 0.5 & & 33 & 9 & 8 & & 0.002 & ** &  &  \\
& 0.55 & & 35 & 6 & 9 & & 1e-05 & **** & MC &  \\
& 0.6 & & 32 & 13 & 5 & & 5e-04 & ** &  &  & (dup) \\
& 0.65 & & 40 & 7 & 3 & & 2e-07 & ***** & MC &  \\
& 0.7 & & 32 & 13 & 5 & & 9e-04 & ** &  &  \\
& 0.75 & & 36 & 13 & 1 & & 5e-04 & ** &  &  \\
\bottomrule
\end{tabular}
\end{table}

