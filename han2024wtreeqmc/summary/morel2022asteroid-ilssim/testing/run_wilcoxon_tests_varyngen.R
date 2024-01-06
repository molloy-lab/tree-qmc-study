#!/usr/bin/env Rscript

library(exactRankTests)
library(coin)

data <- read.csv("../csvs/data-varyngen-for-testing.csv")
data$NGEN <- as.factor(data$NGEN)
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
        for (ngen in ngens) {
        
        xdf <- data[data$NGEN == ngen, ]

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
                         as.character(ngen), "& &",
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
FNR & TREE-QMC-n2 vs. ASTEROID & & 250 & & 12 & 36 & 2 & & 3e-06 & ***** & MC & ( ASTEROID better) \\
FNR & TREE-QMC-n2 vs. ASTEROID & & 500 & & 15 & 28 & 7 & & 0.004 & ** &  & ( ASTEROID better) \\
FNR & TREE-QMC-n2 vs. ASTEROID & & 1000 & & 12 & 27 & 11 & & 0.01 & * &  & ( ASTEROID better) \\
FNR & TREE-QMC-n2 vs. ASTEROID & & 2000 & & 13 & 15 & 22 & & 0.8 &  &  &  \\

FPR & TREE-QMC-n2 vs. ASTEROID & & 250 & & 28 & 22 & 0 & & 0.9 &  &  &  \\
FPR & TREE-QMC-n2 vs. ASTEROID & & 500 & & 26 & 23 & 1 & & 0.3 &  &  &  \\
FPR & TREE-QMC-n2 vs. ASTEROID & & 1000 & & 24 & 23 & 3 & & 0.7 &  &  &  \\
FPR & TREE-QMC-n2 vs. ASTEROID & & 2000 & & 21 & 13 & 16 & & 0.09 &  &  &  \\


FNR & TREE-QMC-n2 vs. ASTER & & 250 & & 28 & 16 & 6 & & 0.1 &  &  &  \\
FNR & TREE-QMC-n2 vs. ASTER & & 500 & & 30 & 14 & 6 & & 0.006 & * &  &  \\
FNR & TREE-QMC-n2 vs. ASTER & & 1000 & & 24 & 19 & 7 & & 0.3 &  &  &  \\
FNR & TREE-QMC-n2 vs. ASTER & & 2000 & & 34 & 9 & 7 & & 0.001 & ** &  &  \\

FPR & TREE-QMC-n2 vs. ASTER & & 250 & & 40 & 10 & 0 & & 2e-08 & ***** & MC &  \\
FPR & TREE-QMC-n2 vs. ASTER & & 500 & & 39 & 11 & 0 & & 3e-07 & ***** & MC &  \\
FPR & TREE-QMC-n2 vs. ASTER & & 1000 & & 28 & 18 & 4 & & 0.009 & * &  &  \\
FPR & TREE-QMC-n2 vs. ASTER & & 2000 & & 39 & 8 & 3 & & 2e-06 & ***** & MC &  \\


FNR & TREE-QMC-n2 vs. wASTRID & & 250 & & 38 & 9 & 3 & & 1e-05 & **** & MC &  \\
FNR & TREE-QMC-n2 vs. wASTRID & & 500 & & 40 & 7 & 3 & & 6e-08 & ***** & MC &  \\
FNR & TREE-QMC-n2 vs. wASTRID & & 1000 & & 35 & 10 & 5 & & 1e-06 & ***** & MC &  \\
FNR & TREE-QMC-n2 vs. wASTRID & & 2000 & & 36 & 8 & 6 & & 3e-07 & ***** & MC &  \\

FPR & TREE-QMC-n2 vs. wASTRID & & 250 & & 45 & 5 & 0 & & 1e-10 & ***** & MC &  \\
FPR & TREE-QMC-n2 vs. wASTRID & & 500 & & 43 & 7 & 0 & & 1e-10 & ***** & MC &  \\
FPR & TREE-QMC-n2 vs. wASTRID & & 1000 & & 39 & 9 & 2 & & 1e-08 & ***** & MC &  \\
FPR & TREE-QMC-n2 vs. wASTRID & & 2000 & & 40 & 6 & 4 & & 1e-09 & ***** & MC &  \\


FNR & TREE-QMC-n2 vs. TREE-QMC-n2-shared & & 250 & & 41 & 5 & 4 & & 4e-10 & ***** & MC &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n2-shared & & 500 & & 41 & 6 & 3 & & 3e-10 & ***** & MC &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n2-shared & & 1000 & & 37 & 7 & 6 & & 5e-08 & ***** & MC &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n2-shared & & 2000 & & 41 & 4 & 5 & & 8e-11 & ***** & MC &  \\

FPR & TREE-QMC-n2 vs. TREE-QMC-n2-shared & & 250 & & 41 & 8 & 1 & & 8e-10 & ***** & MC &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n2-shared & & 500 & & 45 & 5 & 0 & & 9e-11 & ***** & MC &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n2-shared & & 1000 & & 39 & 7 & 4 & & 1e-08 & ***** & MC &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n2-shared & & 2000 & & 42 & 4 & 4 & & 3e-10 & ***** & MC &  \\


FNR & TREE-QMC-n2 vs. TREE-QMC-n1 & & 250 & & 18 & 14 & 18 & & 0.5 &  &  &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n1 & & 500 & & 20 & 13 & 17 & & 0.2 &  &  &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n1 & & 1000 & & 11 & 13 & 26 & & 0.9 &  &  &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n1 & & 2000 & & 14 & 8 & 28 & & 0.2 &  &  &  \\

FPR & TREE-QMC-n2 vs. TREE-QMC-n1 & & 250 & & 16 & 19 & 15 & & 1 &  &  &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n1 & & 500 & & 22 & 16 & 12 & & 0.2 &  &  &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n1 & & 1000 & & 13 & 13 & 24 & & 1 &  &  &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n1 & & 2000 & & 14 & 8 & 28 & & 0.2 &  &  &  \\


FNR & TREE-QMC-n2 vs. TREE-QMC-n0 & & 250 & & 34 & 8 & 8 & & 2e-06 & ***** & MC &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n0 & & 500 & & 33 & 10 & 7 & & 1e-05 & **** & MC &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n0 & & 1000 & & 30 & 12 & 8 & & 6e-04 & ** &  &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n0 & & 2000 & & 36 & 8 & 6 & & 2e-05 & **** & MC &  \\

FPR & TREE-QMC-n2 vs. TREE-QMC-n0 & & 250 & & 35 & 11 & 4 & & 3e-06 & ***** & MC &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n0 & & 500 & & 36 & 9 & 5 & & 4e-06 & ***** & MC &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n0 & & 1000 & & 32 & 13 & 5 & & 5e-04 & ** &  &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n0 & & 2000 & & 36 & 8 & 6 & & 1e-05 & **** & MC &  \\


[1] "done"