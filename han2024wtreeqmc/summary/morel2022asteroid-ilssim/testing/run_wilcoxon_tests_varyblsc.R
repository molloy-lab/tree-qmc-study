#!/usr/bin/env Rscript

library(exactRankTests)
library(coin)

data <- read.csv("../csvs/data-varyblsc-for-testing.csv")
data$BLSC <- as.factor(data$BLSC)
data$REPL <- as.factor(data$REPL)

mthds <- c("ASTEROID", "ASTER", "wASTRID",
           "TREE-QMC-n2-shared", "TREE-QMC-n1", "TREE-QMC-n0")
emets <- c("FNR", "FPR")

psizs <- c(10, 50000000, 100000000, 500000000, 1000000000)  # BASE CONDITION: 50000000
ntaxs <- c(25, 75, 50, 100, 125, 150)  # dup 50
ngens <- c(250, 500, 1000, 2000)  # dup 1000
nbpss <- c(50, 100, 200, 500)  # dup 100
blscs <- c(0.05, 0.1, 1, 10, 100, 200)  # dup 1.0
misss <- c(0.5, 0.55, 0.6, 0.65, 0.7, 0.75)  # dup 0.6

base <- length(psizs)
base <- base + length(ntaxs) - 1
base <- base + length(ngens) - 1
base <- base + length(nbpss) - 1
base <- base + length(blscs) - 1
base <- base + length(misss) - 1
print(paste("There are", base, "model conditions total."))

threshold0 <- 0.05
ntests <- base * length(mthds) * length(emets)
threshold_bonferoni <- threshold0 / ntests

writeme <- paste("alpha =", as.character(threshold0))
print(writeme)
print(paste("bonferroni alpha =", as.character(threshold_bonferoni), 
            "=", as.character(threshold0), "/", as.character(ntests)))

for (mthd in mthds) {
    for (emet in emets) {
        for (blsc in blscs) {
        xdf <- data[data$BLSC == blsc, ]

        # Find ties
        if (emet == "FNR") {
            mthd1 <- xdf$TQMCn2xFNR
            if (mthd == "ASTER") {
                mthd2 <- xdf$ASTERxFNR
            } else if (mthd == "wASTRID") {
                mthd2 <- xdf$WASTRIDxFNR
            } else if (mthd == "ASTEROID") {
                mthd2 <- xdf$ASTEROIDxFNR
            } else if (mthd == "TREE-QMC-n2-shared") {
                mthd2 <- xdf$TQMCn2sharedxFNR
            } else if (mthd == "TREE-QMC-n1") {
                mthd2 <- xdf$TQMCn1xFNR
            } else if (mthd == "TREE-QMC-n0") {
                mthd2 <- xdf$TQMCn0xFNR
            } else {
                print("ERROR")
                exit()
            }
        } else {
            mthd1 <- xdf$TQMCn2xFPR
            if (mthd == "ASTER") {
                mthd2 <- xdf$ASTERxFPR
            } else if (mthd == "wASTRID") {
                mthd2 <- xdf$WASTRIDxFPR
            } else if (mthd == "ASTEROID") {
                mthd2 <- xdf$ASTEROIDxFPR
            } else if (mthd == "TREE-QMC-n2-shared") {
                mthd2 <- xdf$TQMCn2sharedxFPR
            } else if (mthd == "TREE-QMC-n1") {
                mthd2 <- xdf$TQMCn1xFPR
            } else if (mthd == "TREE-QMC-n0") {
                mthd2 <- xdf$TQMCn0xFPR
            } else {
                print("ERROR")
                exit()
            }
        }

        # Count ties
        diff <- mthd2 - mthd1  # positive meants mthd1 is better
        ntie <- sum(diff == 0)
        nother <- sum(diff < 0)
        ntreeqmc <- sum(diff > 0)

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
        note <- ""
        if ((mean(diff) < 0) & (wsr < threshold0)) {
            note <- paste("(", mthd, "better)")
        }

        writeme <- paste(emet, "&", "TREE-QMC-n2 vs.", mthd, 
                        "& &",
                         as.character(blsc), "& &",
                         as.character(ntreeqmc), "&", #"/",
                         as.character(nother), "&", #"/",
                         as.character(ntie), "& &",
                         #format(xwsr, scientific = TRUE), "&", 
                         #format(ywsr, scientific = TRUE), "&",
                         format(wsr, scientific = TRUE, digits=1), "&",
                         stars, "&", mc, "&", note, "\\\\")
        message(writeme)
        }
        message("")
    }
    message("")
}
print("done")
exit()

[1] "There are 26 model conditions total."
[1] "alpha = 0.05"
[1] "bonferroni alpha = 0.00016025641025641 = 0.05 / 312"
FNR & TREE-QMC-n2 vs. ASTEROID & & 0.05 & & 14 & 26 & 10 & & 0.07 &  &  &  \\
FNR & TREE-QMC-n2 vs. ASTEROID & & 0.1 & & 9 & 34 & 7 & & 7e-04 & ** &  & ( ASTEROID better) \\
FNR & TREE-QMC-n2 vs. ASTEROID & & 1 & & 12 & 27 & 11 & & 0.01 & * &  & ( ASTEROID better) \\
FNR & TREE-QMC-n2 vs. ASTEROID & & 10 & & 18 & 18 & 14 & & 0.7 &  &  &  \\
FNR & TREE-QMC-n2 vs. ASTEROID & & 100 & & 21 & 25 & 4 & & 0.3 &  &  &  \\
FNR & TREE-QMC-n2 vs. ASTEROID & & 200 & & 12 & 29 & 9 & & 0.002 & ** &  & ( ASTEROID better) \\

FPR & TREE-QMC-n2 vs. ASTEROID & & 0.05 & & 26 & 23 & 1 & & 0.9 &  &  &  \\
FPR & TREE-QMC-n2 vs. ASTEROID & & 0.1 & & 16 & 29 & 5 & & 0.04 & * &  & ( ASTEROID better) \\
FPR & TREE-QMC-n2 vs. ASTEROID & & 1 & & 24 & 23 & 3 & & 0.7 &  &  &  \\
FPR & TREE-QMC-n2 vs. ASTEROID & & 10 & & 31 & 16 & 3 & & 0.02 & * &  &  \\
FPR & TREE-QMC-n2 vs. ASTEROID & & 100 & & 25 & 23 & 2 & & 0.7 &  &  &  \\
FPR & TREE-QMC-n2 vs. ASTEROID & & 200 & & 17 & 29 & 4 & & 0.02 & * &  & ( ASTEROID better) \\


FNR & TREE-QMC-n2 vs. ASTER & & 0.05 & & 27 & 17 & 6 & & 0.09 &  &  &  \\
FNR & TREE-QMC-n2 vs. ASTER & & 0.1 & & 27 & 15 & 8 & & 0.05 &  &  &  \\
FNR & TREE-QMC-n2 vs. ASTER & & 1 & & 24 & 19 & 7 & & 0.3 &  &  &  \\
FNR & TREE-QMC-n2 vs. ASTER & & 10 & & 29 & 13 & 8 & & 6e-04 & ** &  &  \\
FNR & TREE-QMC-n2 vs. ASTER & & 100 & & 37 & 10 & 3 & & 2e-06 & ***** & MC &  \\
FNR & TREE-QMC-n2 vs. ASTER & & 200 & & 25 & 21 & 4 & & 0.6 &  &  &  \\

FPR & TREE-QMC-n2 vs. ASTER & & 0.05 & & 32 & 16 & 2 & & 0.002 & ** &  &  \\
FPR & TREE-QMC-n2 vs. ASTER & & 0.1 & & 29 & 15 & 6 & & 0.002 & ** &  &  \\
FPR & TREE-QMC-n2 vs. ASTER & & 1 & & 28 & 18 & 4 & & 0.009 & * &  &  \\
FPR & TREE-QMC-n2 vs. ASTER & & 10 & & 37 & 12 & 1 & & 5e-07 & ***** & MC &  \\
FPR & TREE-QMC-n2 vs. ASTER & & 100 & & 41 & 9 & 0 & & 8e-08 & ***** & MC &  \\
FPR & TREE-QMC-n2 vs. ASTER & & 200 & & 28 & 20 & 2 & & 0.2 &  &  &  \\


FNR & TREE-QMC-n2 vs. wASTRID & & 0.05 & & 30 & 14 & 6 & & 0.005 & * &  &  \\
FNR & TREE-QMC-n2 vs. wASTRID & & 0.1 & & 35 & 10 & 5 & & 1e-05 & **** & MC &  \\
FNR & TREE-QMC-n2 vs. wASTRID & & 1 & & 35 & 10 & 5 & & 1e-06 & ***** & MC &  \\
FNR & TREE-QMC-n2 vs. wASTRID & & 10 & & 34 & 13 & 3 & & 3e-05 & **** & MC &  \\
FNR & TREE-QMC-n2 vs. wASTRID & & 100 & & 34 & 14 & 2 & & 7e-04 & ** &  &  \\
FNR & TREE-QMC-n2 vs. wASTRID & & 200 & & 26 & 23 & 1 & & 0.6 &  &  &  \\

FPR & TREE-QMC-n2 vs. wASTRID & & 0.05 & & 35 & 14 & 1 & & 3e-04 & *** &  &  \\
FPR & TREE-QMC-n2 vs. wASTRID & & 0.1 & & 38 & 10 & 2 & & 2e-07 & ***** & MC &  \\
FPR & TREE-QMC-n2 vs. wASTRID & & 1 & & 39 & 9 & 2 & & 1e-08 & ***** & MC &  \\
FPR & TREE-QMC-n2 vs. wASTRID & & 10 & & 39 & 10 & 1 & & 1e-07 & ***** & MC &  \\
FPR & TREE-QMC-n2 vs. wASTRID & & 100 & & 35 & 14 & 1 & & 2e-04 & *** &  &  \\
FPR & TREE-QMC-n2 vs. wASTRID & & 200 & & 27 & 23 & 0 & & 0.3 &  &  &  \\


FNR & TREE-QMC-n2 vs. TREE-QMC-n2-shared & & 0.05 & & 49 & 1 & 0 & & 5e-15 & ***** & MC &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n2-shared & & 0.1 & & 48 & 2 & 0 & & 1e-13 & ***** & MC &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n2-shared & & 1 & & 37 & 7 & 6 & & 5e-08 & ***** & MC &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n2-shared & & 10 & & 43 & 6 & 1 & & 3e-10 & ***** & MC &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n2-shared & & 100 & & 48 & 2 & 0 & & 2e-14 & ***** & MC &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n2-shared & & 200 & & 46 & 4 & 0 & & 1e-11 & ***** & MC &  \\

FPR & TREE-QMC-n2 vs. TREE-QMC-n2-shared & & 0.05 & & 49 & 1 & 0 & & 4e-15 & ***** & MC &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n2-shared & & 0.1 & & 48 & 2 & 0 & & 6e-14 & ***** & MC &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n2-shared & & 1 & & 39 & 7 & 4 & & 1e-08 & ***** & MC &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n2-shared & & 10 & & 44 & 6 & 0 & & 6e-11 & ***** & MC &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n2-shared & & 100 & & 49 & 1 & 0 & & 1e-14 & ***** & MC &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n2-shared & & 200 & & 47 & 3 & 0 & & 9e-12 & ***** & MC &  \\


FNR & TREE-QMC-n2 vs. TREE-QMC-n1 & & 0.05 & & 15 & 24 & 11 & & 0.3 &  &  &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n1 & & 0.1 & & 14 & 22 & 14 & & 0.2 &  &  &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n1 & & 1 & & 11 & 13 & 26 & & 0.9 &  &  &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n1 & & 10 & & 21 & 10 & 19 & & 0.01 & * &  &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n1 & & 100 & & 30 & 12 & 8 & & 0.002 & ** &  &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n1 & & 200 & & 22 & 19 & 9 & & 0.6 &  &  &  \\

FPR & TREE-QMC-n2 vs. TREE-QMC-n1 & & 0.05 & & 16 & 25 & 9 & & 0.4 &  &  &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n1 & & 0.1 & & 14 & 23 & 13 & & 0.2 &  &  &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n1 & & 1 & & 13 & 13 & 24 & & 1 &  &  &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n1 & & 10 & & 24 & 10 & 16 & & 0.02 & * &  &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n1 & & 100 & & 31 & 12 & 7 & & 0.01 & * &  &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n1 & & 200 & & 22 & 19 & 9 & & 0.5 &  &  &  \\


FNR & TREE-QMC-n2 vs. TREE-QMC-n0 & & 0.05 & & 29 & 14 & 7 & & 0.01 & * &  &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n0 & & 0.1 & & 34 & 10 & 6 & & 7e-05 & *** & MC &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n0 & & 1 & & 30 & 12 & 8 & & 6e-04 & ** &  &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n0 & & 10 & & 34 & 12 & 4 & & 1e-05 & **** & MC &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n0 & & 100 & & 37 & 7 & 6 & & 3e-08 & ***** & MC &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n0 & & 200 & & 31 & 13 & 6 & & 6e-04 & ** &  &  \\

FPR & TREE-QMC-n2 vs. TREE-QMC-n0 & & 0.05 & & 33 & 15 & 2 & & 0.004 & ** &  &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n0 & & 0.1 & & 35 & 12 & 3 & & 3e-05 & **** & MC &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n0 & & 1 & & 32 & 13 & 5 & & 5e-04 & ** &  &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n0 & & 10 & & 37 & 12 & 1 & & 5e-06 & ***** & MC &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n0 & & 100 & & 40 & 8 & 2 & & 3e-08 & ***** & MC &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n0 & & 200 & & 35 & 13 & 2 & & 3e-04 & *** &  &  \\


[1] "done"


[1] "done"