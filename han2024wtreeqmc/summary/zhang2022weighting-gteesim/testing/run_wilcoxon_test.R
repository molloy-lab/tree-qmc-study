#!/usr/bin/env Rscript

# install.packages("coin", dependencies = TRUE)

library(exactRankTests)
library(coin)

#https://sphweb.bumc.bu.edu/otlt/mph-modules/bs/bs704_nonparametric/BS704_Nonparametric6.html
#https://stat.ethz.ch/pipermail/r-help/2011-April/274931.html
#https://www.graphpad.com/guides/prism/latest/statistics/stat_pratt.htm
#https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/wilcox.test
# https://www.rdocumentation.org/packages/rstatix/versions/0.7.2/topics/wilcox_test

data <- read.csv("../csvs/data-for-testing.csv")
data$NBPS <- as.factor(data$NBPS)
data$NGEN <- as.factor(data$NGEN)
data$SUPP <- as.factor(data$SUPP)
data$REPL <- as.factor(data$REPL)

ngens <- c("50", "200", "500", "1000")
nbpss <- c("200", "400", "800", "1600")
mthds <- c("ASTER-wh", "ASTRID-ws", "TREE-QMC-n2")
supps <- c("abayes", "bs")

threshold0 <- 0.05
ntests <- length(ngens) * length(nbpss) * length(mthds) * length(supps)
threshold_bonferoni <- threshold0 / ntests

#writeme <- paste("alpha =", as.character(threshold0))
#print(writeme)
#print(paste("bonferroni alpha =", as.character(threshold_bonferoni), 
#            "=", as.character(threshold0), "/", as.character(ntests)))

for (mthd in mthds) {

    message("\\begin{table}[!h]")
    message(paste("\\caption[Statistical testing for TREE-QMC-wh\\_n2 vs", mthd)) 
    message("on S100 data]{\\textbf{Testing for differences between TREE-QMC-wh\\_n2 vs")
    message(paste(mthd, "simulated data sets.} BET is the number of replicates for which"))
    message(paste("TREE-QMC-wh\\_n2 achieves a lower RF error rate and thus is better than", mthd)) 
    message("WOR is the number of replicates for which wTREE-QMC-wh\\_n2  achieves a higher")
    message(paste("RF error rate and thus is worse than", mthd, ",")) 
    message("and TIE is the number of replicates where the two methods tie.")
    message("Significance is evaluated using Wilcoxon signed-rank tests on the species tree (RF) error rates.")
    message("The symbols *, **, ***, ****, ***** indicate significance at $p <$ 0.5, 0.005,$ and so on.")
    message(paste("Note: $p < ", as.character(threshold0), "/", as.character(ntests), "=", format(threshold_bonferoni, scientific = TRUE, digits=1)))
    message(paste("is be significant after Bonferroni correction for the ", as.character(ntests), " tests made on the S100 data sets (MC).}"))
    if (mthd == "TREE-QMC-n2") {
        message("IMPORTANT: All other results for TREE-QMC-n2 do NOT use the refined tree from IQTREE.")
    }
    message("\\centering")
    message("\\normal")
    message("\\begin{tabular}{r r l r r r l c c c l}")
    message("\\toprule") 
    message(paste("\\multicolumn{11}{c}{\\textbf{TREE-QMC-wh\\_n2  vs ", mthd, "}} \\\\"))
    message("\\midrule")
    message("\\# of genes & sequence length & & BET & WOR & TIE & & p-val & sig & MC & note \\\\")

    for (supp in supps) {
        message("\\midrule")
        if (supp == "bs") {
            if (mthd == "TREE-QMC-n2") { okgo <- "TREE-QMC used to randomly refined polytomies" } 
            else { okgo <- "Boostrap Support" }
        } else {
            if (mthd == "TREE-QMC-n2") { okgo <- "IQtree used refined polytomies (part of abayes)" } 
            else { okgo <- "Abayes Support" }
        }
        message(paste("\\multicolumn{11}{c}{\\textbf{", okgo, "}} \\\\"))
        message("\\midrule")

        for (ngen in ngens) {
            for (nbps in nbpss) {
        
        df <- data[data$SUPP == supp, ]
        df <- df[df$NGEN == ngen, ]
        df <- df[df$NBPS == nbps, ]

        # Set-up
        mthd1 <- df$TQMCwhn2xSERF #FN
        if (mthd == "ASTER-wh") {
            mthd2 <- df$ASTERHxSERF #FN
            supplab = supp
        } else if (mthd == "ASTRID-ws") {
            mthd2 <- df$WASTRIDxSERF #FN
            supplab = supp
        } else if (mthd == "TREE-QMC-n2") {
            mthd2 <- df$TQMCn2xSERF #FN
            #print(mean(mthd2))
            supplab = supp #"none"
        } else {
            print("ERROR")
            exit()
        }

        # Count ties
        diff <- mthd2 - mthd1  # lower is better
        ntie <- sum(diff == 0)
        nother <- sum(diff < 0)
        ntreeqmc <- sum(diff > 0)

        # Perform wilcoxon test
        #xwsr <-  wilcox.exact(mthd1, mthd2,
        #                      paired = TRUE,
        #                      alternative = "two.sided",
        #                      conf.int = TRUE)
        #xwsr <- xwsr$p.value
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
        note <- ""
        if ((mean(diff) < 0) & (wsr < threshold0)) {
            note <- paste("(", mthd, "better)")
        }

        #writeme <- paste("TQMC-wh_n2 vs.", mthd, "&",
        #                 supplab, "& &",
        #                 as.character(ngen), "&",
        #                 as.character(nbps), "& &",
        #                 as.character(ntreeqmc), "&",
        #                 as.character(nother), "&",
        #                 as.character(ntie), "& &",
        #                 format(wsr, scientific = TRUE, digits=1), "&", 
        #                 stars, "&", mc, "&", note, " \\\\")
        writeme <- paste(as.character(ngen), "&",
                         as.character(nbps), "& &",
                         as.character(ntreeqmc), "&",
                         as.character(nother), "&",
                         as.character(ntie), "& &",
                         format(wsr, scientific = TRUE, digits=1), "&", 
                         stars, "&", mc, "&", note, " \\\\")
        message(writeme)

            }
        }
    }
    message("\\bottomrule")
    message("\\end{tabular}")
    message("\\end{table}")
    message("")
}

print("done")
exit()

\begin{table}[!h]
\caption[Statistical testing for TREE-QMC-wh\_n2 vs ASTER-wh
on S100 data]{\textbf{Testing for differences between TREE-QMC-wh\_n2 vs
ASTER-wh simulated data sets.} BET is the number of replicates for which
TREE-QMC-wh\_n2 achieves a lower RF error rate and thus is better than ASTER-wh
WOR is the number of replicates for which wTREE-QMC-wh\_n2  achieves a higher
RF error rate and thus is worse than ASTER-wh ,
and TIE is the number of replicates where the two methods tie.
Significance is evaluated using Wilcoxon signed-rank tests on the species tree (RF) error rates.
The symbols *, **, ***, ****, ***** indicate significance at $p <$ 0.5, 0.005,$ and so on.
Note: $p <  0.05 / 96 = 5e-04
is be significant after Bonferroni correction for the  96  tests made on the S100 data sets (MC).}
\centering
\normal
\begin{tabular}{r r l r r r l c c c l}
\toprule
\multicolumn{11}{c}{\textbf{TREE-QMC-wh\_n2  vs  ASTER-wh }} \\
\midrule
\# of genes & sequence length & & BET & WOR & TIE & & p-val & sig & MC & note \\
\midrule
\multicolumn{11}{c}{\textbf{ Abayes Support }} \\
\midrule
50 & 200 & & 23 & 17 & 10 & & 0.3 &  &  &   \\
50 & 400 & & 33 & 7 & 10 & & 4e-04 & *** & MC &   \\
50 & 800 & & 27 & 8 & 15 & & 3e-04 & *** & MC &   \\
50 & 1600 & & 23 & 6 & 21 & & 4e-04 & *** & MC &   \\
200 & 200 & & 30 & 9 & 11 & & 0.007 & * &  &   \\
200 & 400 & & 21 & 14 & 15 & & 0.2 &  &  &   \\
200 & 800 & & 18 & 11 & 21 & & 0.2 &  &  &   \\
200 & 1600 & & 22 & 10 & 18 & & 0.02 & * &  &   \\
500 & 200 & & 22 & 15 & 13 & & 0.8 &  &  &   \\
500 & 400 & & 14 & 22 & 14 & & 0.5 &  &  &   \\
500 & 800 & & 18 & 11 & 21 & & 0.5 &  &  &   \\
500 & 1600 & & 21 & 13 & 16 & & 0.2 &  &  &   \\
1000 & 200 & & 20 & 9 & 21 & & 0.006 & * &  &   \\
1000 & 400 & & 20 & 6 & 24 & & 0.003 & ** &  &   \\
1000 & 800 & & 20 & 7 & 23 & & 0.008 & * &  &   \\
1000 & 1600 & & 19 & 5 & 26 & & 0.002 & ** &  &   \\
\midrule
\multicolumn{11}{c}{\textbf{ Boostrap Support }} \\
\midrule
50 & 200 & & 22 & 17 & 11 & & 0.4 &  &  &   \\
50 & 400 & & 19 & 14 & 17 & & 0.5 &  &  &   \\
50 & 800 & & 23 & 18 & 9 & & 0.8 &  &  &   \\
50 & 1600 & & 21 & 16 & 13 & & 0.2 &  &  &   \\
200 & 200 & & 19 & 17 & 14 & & 0.6 &  &  &   \\
200 & 400 & & 21 & 11 & 18 & & 0.06 &  &  &   \\
200 & 800 & & 17 & 18 & 15 & & 0.6 &  &  &   \\
200 & 1600 & & 17 & 11 & 22 & & 0.2 &  &  &   \\
500 & 200 & & 23 & 18 & 9 & & 0.8 &  &  &   \\
500 & 400 & & 21 & 14 & 15 & & 0.1 &  &  &   \\
500 & 800 & & 17 & 10 & 23 & & 0.3 &  &  &   \\
500 & 1600 & & 13 & 8 & 29 & & 0.2 &  &  &   \\
1000 & 200 & & 16 & 19 & 15 & & 0.6 &  &  &   \\
1000 & 400 & & 16 & 10 & 24 & & 0.1 &  &  &   \\
1000 & 800 & & 17 & 6 & 27 & & 0.09 &  &  &   \\
1000 & 1600 & & 13 & 10 & 27 & & 0.4 &  &  &   \\
\bottomrule
\end{tabular}
\end{table}

\begin{table}[!h]
\caption[Statistical testing for TREE-QMC-wh\_n2 vs ASTRID-ws
on S100 data]{\textbf{Testing for differences between TREE-QMC-wh\_n2 vs
ASTRID-ws simulated data sets.} BET is the number of replicates for which
TREE-QMC-wh\_n2 achieves a lower RF error rate and thus is better than ASTRID-ws
WOR is the number of replicates for which wTREE-QMC-wh\_n2  achieves a higher
RF error rate and thus is worse than ASTRID-ws ,
and TIE is the number of replicates where the two methods tie.
Significance is evaluated using Wilcoxon signed-rank tests on the species tree (RF) error rates.
The symbols *, **, ***, ****, ***** indicate significance at $p <$ 0.5, 0.005,$ and so on.
Note: $p <  0.05 / 96 = 5e-04
is be significant after Bonferroni correction for the  96  tests made on the S100 data sets (MC).}
\centering
\normal
\begin{tabular}{r r l r r r l c c c l}
\toprule
\multicolumn{11}{c}{\textbf{TREE-QMC-wh\_n2  vs  ASTRID-ws }} \\
\midrule
\# of genes & sequence length & & BET & WOR & TIE & & p-val & sig & MC & note \\
\midrule
\multicolumn{11}{c}{\textbf{ Abayes Support }} \\
\midrule
50 & 200 & & 30 & 16 & 4 & & 0.03 & * &  &   \\
50 & 400 & & 33 & 12 & 5 & & 0.001 & ** &  &   \\
50 & 800 & & 35 & 10 & 5 & & 2e-04 & *** & MC &   \\
50 & 1600 & & 34 & 10 & 6 & & 2e-05 & **** & MC &   \\
200 & 200 & & 27 & 7 & 16 & & 3e-04 & *** & MC &   \\
200 & 400 & & 25 & 13 & 12 & & 0.03 & * &  &   \\
200 & 800 & & 25 & 13 & 12 & & 0.03 & * &  &   \\
200 & 1600 & & 29 & 7 & 14 & & 0.006 & * &  &   \\
500 & 200 & & 24 & 16 & 10 & & 0.1 &  &  &   \\
500 & 400 & & 23 & 13 & 14 & & 0.1 &  &  &   \\
500 & 800 & & 22 & 11 & 17 & & 0.02 & * &  &   \\
500 & 1600 & & 24 & 9 & 17 & & 0.04 & * &  &   \\
1000 & 200 & & 25 & 9 & 16 & & 0.008 & * &  &   \\
1000 & 400 & & 21 & 13 & 16 & & 0.07 &  &  &   \\
1000 & 800 & & 21 & 8 & 21 & & 0.004 & ** &  &   \\
1000 & 1600 & & 23 & 5 & 22 & & 4e-04 & *** & MC &   \\
\midrule
\multicolumn{11}{c}{\textbf{ Boostrap Support }} \\
\midrule
50 & 200 & & 30 & 12 & 8 & & 0.004 & ** &  &   \\
50 & 400 & & 31 & 12 & 7 & & 0.002 & ** &  &   \\
50 & 800 & & 31 & 13 & 6 & & 0.002 & ** &  &   \\
50 & 1600 & & 32 & 9 & 9 & & 2e-04 & *** & MC &   \\
200 & 200 & & 20 & 18 & 12 & & 0.6 &  &  &   \\
200 & 400 & & 22 & 14 & 14 & & 0.03 & * &  &   \\
200 & 800 & & 23 & 11 & 16 & & 0.02 & * &  &   \\
200 & 1600 & & 19 & 11 & 20 & & 0.04 & * &  &   \\
500 & 200 & & 17 & 14 & 19 & & 0.6 &  &  &   \\
500 & 400 & & 25 & 9 & 16 & & 0.002 & ** &  &   \\
500 & 800 & & 17 & 12 & 21 & & 0.2 &  &  &   \\
500 & 1600 & & 20 & 9 & 21 & & 0.03 & * &  &   \\
1000 & 200 & & 23 & 18 & 9 & & 1 &  &  &   \\
1000 & 400 & & 21 & 11 & 18 & & 0.09 &  &  &   \\
1000 & 800 & & 15 & 9 & 26 & & 0.1 &  &  &   \\
1000 & 1600 & & 15 & 12 & 23 & & 0.4 &  &  &   \\
\bottomrule
\end{tabular}
\end{table}

\begin{table}[!h]
\caption[Statistical testing for TREE-QMC-wh\_n2 vs TREE-QMC-n2
on S100 data]{\textbf{Testing for differences between TREE-QMC-wh\_n2 vs
TREE-QMC-n2 simulated data sets.} BET is the number of replicates for which
TREE-QMC-wh\_n2 achieves a lower RF error rate and thus is better than TREE-QMC-n2
WOR is the number of replicates for which wTREE-QMC-wh\_n2  achieves a higher
RF error rate and thus is worse than TREE-QMC-n2 ,
and TIE is the number of replicates where the two methods tie.
Significance is evaluated using Wilcoxon signed-rank tests on the species tree (RF) error rates.
The symbols *, **, ***, ****, ***** indicate significance at $p <$ 0.5, 0.005,$ and so on.
Note: $p <  0.05 / 96 = 5e-04
is be significant after Bonferroni correction for the  96  tests made on the S100 data sets (MC).}
IMPORTANT: All other results for TREE-QMC-n2 do NOT use the refined tree from IQTREE.
\centering
\normal
\begin{tabular}{r r l r r r l c c c l}
\toprule
\multicolumn{11}{c}{\textbf{TREE-QMC-wh\_n2  vs  TREE-QMC-n2 }} \\
\midrule
\# of genes & sequence length & & BET & WOR & TIE & & p-val & sig & MC & note \\
\midrule
\multicolumn{11}{c}{\textbf{ IQtree used refined polytomies (part of abayes) }} \\
\midrule
50 & 200 & & 38 & 8 & 4 & & 3e-07 & ***** & MC &   \\
50 & 400 & & 31 & 11 & 8 & & 1e-04 & *** & MC &   \\
50 & 800 & & 32 & 8 & 10 & & 7e-05 & *** & MC &   \\
50 & 1600 & & 24 & 16 & 10 & & 0.05 &  &  &   \\
200 & 200 & & 35 & 11 & 4 & & 1e-06 & ***** & MC &   \\
200 & 400 & & 34 & 8 & 8 & & 7e-06 & **** & MC &   \\
200 & 800 & & 22 & 11 & 17 & & 0.03 & * &  &   \\
200 & 1600 & & 29 & 9 & 12 & & 9e-04 & ** &  &   \\
500 & 200 & & 29 & 9 & 12 & & 6e-04 & ** &  &   \\
500 & 400 & & 28 & 9 & 13 & & 1e-04 & *** & MC &   \\
500 & 800 & & 22 & 9 & 19 & & 0.008 & * &  &   \\
500 & 1600 & & 24 & 6 & 20 & & 2e-04 & *** & MC &   \\
1000 & 200 & & 26 & 6 & 18 & & 4e-05 & **** & MC &   \\
1000 & 400 & & 29 & 4 & 17 & & 8e-07 & ***** & MC &   \\
1000 & 800 & & 28 & 5 & 17 & & 9e-06 & **** & MC &   \\
1000 & 1600 & & 22 & 5 & 23 & & 0.002 & ** &  &   \\
\midrule
\multicolumn{11}{c}{\textbf{ TREE-QMC used to randomly refined polytomies }} \\
\midrule
50 & 200 & & 29 & 13 & 8 & & 3e-04 & *** & MC &   \\
50 & 400 & & 25 & 20 & 5 & & 0.4 &  &  &   \\
50 & 800 & & 26 & 15 & 9 & & 0.04 & * &  &   \\
50 & 1600 & & 22 & 18 & 10 & & 0.3 &  &  &   \\
200 & 200 & & 24 & 16 & 10 & & 0.05 & * &  &   \\
200 & 400 & & 29 & 12 & 9 & & 0.005 & * &  &   \\
200 & 800 & & 17 & 14 & 19 & & 0.3 &  &  &   \\
200 & 1600 & & 24 & 9 & 17 & & 0.003 & ** &  &   \\
500 & 200 & & 24 & 19 & 7 & & 0.05 & * &  &   \\
500 & 400 & & 26 & 9 & 15 & & 2e-04 & *** & MC &   \\
500 & 800 & & 18 & 16 & 16 & & 0.3 &  &  &   \\
500 & 1600 & & 17 & 12 & 21 & & 0.2 &  &  &   \\
1000 & 200 & & 26 & 17 & 7 & & 0.04 & * &  &   \\
1000 & 400 & & 27 & 10 & 13 & & 0.001 & ** &  &   \\
1000 & 800 & & 21 & 9 & 20 & & 0.005 & ** &  &   \\
1000 & 1600 & & 12 & 11 & 27 & & 0.6 &  &  &   \\
\bottomrule
\end{tabular}
\end{table}


