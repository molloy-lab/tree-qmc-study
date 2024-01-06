#!/usr/bin/env Rscript

library(exactRankTests)
library(coin)

data <- read.csv("../csvs/data-varynbps-for-testing.csv")
data$NBPS <- as.factor(data$NBPS)
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

threshold0 <- 0.05
ntests <- base * length(mthds) * length(emets)
threshold_bonferoni <- threshold0 / ntests

writeme <- paste("alpha =", as.character(threshold0))
print(writeme)
print(paste("bonferroni alpha =", as.character(threshold_bonferoni), 
            "=", as.character(threshold0), "/", as.character(ntests)))

for (mthd in mthds) {
    for (emet in emets) {
        for (nbps in nbpss) {
        
        xdf <- data[data$NBPS == nbps, ]

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
                         as.character(nbps), "& &",
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

[1] "alpha = 0.05"
[1] "bonferroni alpha = 0.00016025641025641 = 0.05 / 312"
FNR & TREE-QMC-n2 vs. ASTEROID & & 50 & & 13 & 28 & 9 & & 0.02 & * &  & ( ASTEROID better) \\
FNR & TREE-QMC-n2 vs. ASTEROID & & 100 & & 12 & 27 & 11 & & 0.01 & * &  & ( ASTEROID better) \\
FNR & TREE-QMC-n2 vs. ASTEROID & & 200 & & 13 & 24 & 13 & & 0.03 & * &  & ( ASTEROID better) \\
FNR & TREE-QMC-n2 vs. ASTEROID & & 500 & & 15 & 20 & 15 & & 0.3 &  &  &  \\

FPR & TREE-QMC-n2 vs. ASTEROID & & 50 & & 19 & 27 & 4 & & 0.3 &  &  &  \\
FPR & TREE-QMC-n2 vs. ASTEROID & & 100 & & 24 & 23 & 3 & & 0.7 &  &  &  \\
FPR & TREE-QMC-n2 vs. ASTEROID & & 200 & & 27 & 19 & 4 & & 0.2 &  &  &  \\
FPR & TREE-QMC-n2 vs. ASTEROID & & 500 & & 32 & 11 & 7 & & 8e-04 & ** &  &  \\


FNR & TREE-QMC-n2 vs. ASTER & & 50 & & 29 & 11 & 10 & & 0.001 & ** &  &  \\
FNR & TREE-QMC-n2 vs. ASTER & & 100 & & 24 & 19 & 7 & & 0.3 &  &  &  \\
FNR & TREE-QMC-n2 vs. ASTER & & 200 & & 17 & 26 & 7 & & 0.4 &  &  &  \\
FNR & TREE-QMC-n2 vs. ASTER & & 500 & & 21 & 17 & 12 & & 0.3 &  &  &  \\

FPR & TREE-QMC-n2 vs. ASTER & & 50 & & 36 & 10 & 4 & & 2e-05 & **** & MC &  \\
FPR & TREE-QMC-n2 vs. ASTER & & 100 & & 28 & 18 & 4 & & 0.009 & * &  &  \\
FPR & TREE-QMC-n2 vs. ASTER & & 200 & & 29 & 18 & 3 & & 0.009 & * &  &  \\
FPR & TREE-QMC-n2 vs. ASTER & & 500 & & 33 & 14 & 3 & & 2e-04 & *** &  &  \\


FNR & TREE-QMC-n2 vs. wASTRID & & 50 & & 30 & 17 & 3 & & 0.007 & * &  &  \\
FNR & TREE-QMC-n2 vs. wASTRID & & 100 & & 35 & 10 & 5 & & 1e-06 & ***** & MC &  \\
FNR & TREE-QMC-n2 vs. wASTRID & & 200 & & 41 & 6 & 3 & & 9e-08 & ***** & MC &  \\
FNR & TREE-QMC-n2 vs. wASTRID & & 500 & & 42 & 4 & 4 & & 6e-10 & ***** & MC &  \\

FPR & TREE-QMC-n2 vs. wASTRID & & 50 & & 31 & 17 & 2 & & 3e-04 & *** &  &  \\
FPR & TREE-QMC-n2 vs. wASTRID & & 100 & & 39 & 9 & 2 & & 1e-08 & ***** & MC &  \\
FPR & TREE-QMC-n2 vs. wASTRID & & 200 & & 45 & 5 & 0 & & 1e-10 & ***** & MC &  \\
FPR & TREE-QMC-n2 vs. wASTRID & & 500 & & 45 & 2 & 3 & & 4e-12 & ***** & MC &  \\


FNR & TREE-QMC-n2 vs. TREE-QMC-n2-shared & & 50 & & 41 & 3 & 6 & & 2e-11 & ***** & MC &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n2-shared & & 100 & & 37 & 7 & 6 & & 5e-08 & ***** & MC &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n2-shared & & 200 & & 37 & 7 & 6 & & 9e-08 & ***** & MC &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n2-shared & & 500 & & 40 & 4 & 6 & & 4e-08 & ***** & MC &  \\

FPR & TREE-QMC-n2 vs. TREE-QMC-n2-shared & & 50 & & 42 & 6 & 2 & & 1e-11 & ***** & MC &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n2-shared & & 100 & & 39 & 7 & 4 & & 1e-08 & ***** & MC &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n2-shared & & 200 & & 38 & 7 & 5 & & 1e-07 & ***** & MC &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n2-shared & & 500 & & 40 & 6 & 4 & & 1e-08 & ***** & MC &  \\


FNR & TREE-QMC-n2 vs. TREE-QMC-n1 & & 50 & & 16 & 17 & 17 & & 0.5 &  &  &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n1 & & 100 & & 11 & 13 & 26 & & 0.9 &  &  &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n1 & & 200 & & 14 & 14 & 22 & & 1 &  &  &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n1 & & 500 & & 6 & 13 & 31 & & 0.09 &  &  &  \\

FPR & TREE-QMC-n2 vs. TREE-QMC-n1 & & 50 & & 16 & 17 & 17 & & 0.7 &  &  &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n1 & & 100 & & 13 & 13 & 24 & & 1 &  &  &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n1 & & 200 & & 14 & 14 & 22 & & 0.8 &  &  &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n1 & & 500 & & 6 & 13 & 31 & & 0.06 &  &  &  \\


FNR & TREE-QMC-n2 vs. TREE-QMC-n0 & & 50 & & 31 & 9 & 10 & & 1e-04 & *** & MC &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n0 & & 100 & & 30 & 12 & 8 & & 6e-04 & ** &  &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n0 & & 200 & & 26 & 11 & 13 & & 0.004 & ** &  &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n0 & & 500 & & 27 & 8 & 15 & & 0.002 & ** &  &  \\

FPR & TREE-QMC-n2 vs. TREE-QMC-n0 & & 50 & & 32 & 9 & 9 & & 5e-05 & **** & MC &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n0 & & 100 & & 32 & 13 & 5 & & 5e-04 & ** &  &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n0 & & 200 & & 28 & 13 & 9 & & 0.004 & ** &  &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n0 & & 500 & & 28 & 9 & 13 & & 4e-04 & *** &  &  \\


[1] "done"