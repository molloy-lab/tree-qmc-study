#!/usr/bin/env Rscript

library(exactRankTests)
library(coin)

data <- read.csv("../csvs/data-varymiss-for-testing.csv")
data$MISS <- as.factor(data$MISS)
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
        for (miss in misss) {
        
        xdf <- data[data$MISS == miss, ]

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
                         as.character(miss), "& &",
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
FNR & TREE-QMC-n2 vs. ASTEROID & & 0.5 & & 14 & 23 & 13 & & 0.3 &  &  &  \\
FNR & TREE-QMC-n2 vs. ASTEROID & & 0.55 & & 15 & 17 & 18 & & 0.8 &  &  &  \\
FNR & TREE-QMC-n2 vs. ASTEROID & & 0.6 & & 12 & 27 & 11 & & 0.01 & * &  & ( ASTEROID better) \\
FNR & TREE-QMC-n2 vs. ASTEROID & & 0.65 & & 15 & 24 & 11 & & 0.08 &  &  &  \\
FNR & TREE-QMC-n2 vs. ASTEROID & & 0.7 & & 10 & 32 & 8 & & 5e-06 & ***** & MC & ( ASTEROID better) \\
FNR & TREE-QMC-n2 vs. ASTEROID & & 0.75 & & 7 & 40 & 3 & & 5e-09 & ***** & MC & ( ASTEROID better) \\

FPR & TREE-QMC-n2 vs. ASTEROID & & 0.5 & & 17 & 23 & 10 & & 0.3 &  &  &  \\
FPR & TREE-QMC-n2 vs. ASTEROID & & 0.55 & & 22 & 16 & 12 & & 0.3 &  &  &  \\
FPR & TREE-QMC-n2 vs. ASTEROID & & 0.6 & & 24 & 23 & 3 & & 0.7 &  &  &  \\
FPR & TREE-QMC-n2 vs. ASTEROID & & 0.65 & & 39 & 11 & 0 & & 1e-05 & **** & MC &  \\
FPR & TREE-QMC-n2 vs. ASTEROID & & 0.7 & & 33 & 17 & 0 & & 0.07 &  &  &  \\
FPR & TREE-QMC-n2 vs. ASTEROID & & 0.75 & & 28 & 22 & 0 & & 0.2 &  &  &  \\


FNR & TREE-QMC-n2 vs. ASTER & & 0.5 & & 27 & 11 & 12 & & 0.02 & * &  &  \\
FNR & TREE-QMC-n2 vs. ASTER & & 0.55 & & 31 & 9 & 10 & & 5e-04 & *** &  &  \\
FNR & TREE-QMC-n2 vs. ASTER & & 0.6 & & 24 & 19 & 7 & & 0.3 &  &  &  \\
FNR & TREE-QMC-n2 vs. ASTER & & 0.65 & & 29 & 7 & 14 & & 2e-04 & *** & MC &  \\
FNR & TREE-QMC-n2 vs. ASTER & & 0.7 & & 22 & 20 & 8 & & 0.4 &  &  &  \\
FNR & TREE-QMC-n2 vs. ASTER & & 0.75 & & 28 & 15 & 7 & & 0.01 & * &  &  \\

FPR & TREE-QMC-n2 vs. ASTER & & 0.5 & & 28 & 11 & 11 & & 0.009 & * &  &  \\
FPR & TREE-QMC-n2 vs. ASTER & & 0.55 & & 36 & 9 & 5 & & 3e-05 & **** & MC &  \\
FPR & TREE-QMC-n2 vs. ASTER & & 0.6 & & 28 & 18 & 4 & & 0.009 & * &  &  \\
FPR & TREE-QMC-n2 vs. ASTER & & 0.65 & & 45 & 4 & 1 & & 1e-12 & ***** & MC &  \\
FPR & TREE-QMC-n2 vs. ASTER & & 0.7 & & 43 & 7 & 0 & & 8e-09 & ***** & MC &  \\
FPR & TREE-QMC-n2 vs. ASTER & & 0.75 & & 45 & 5 & 0 & & 3e-10 & ***** & MC &  \\


FNR & TREE-QMC-n2 vs. wASTRID & & 0.5 & & 29 & 10 & 11 & & 4e-04 & *** &  &  \\
FNR & TREE-QMC-n2 vs. wASTRID & & 0.55 & & 32 & 15 & 3 & & 0.003 & ** &  &  \\
FNR & TREE-QMC-n2 vs. wASTRID & & 0.6 & & 35 & 10 & 5 & & 1e-06 & ***** & MC &  \\
FNR & TREE-QMC-n2 vs. wASTRID & & 0.65 & & 35 & 9 & 6 & & 4e-08 & ***** & MC &  \\
FNR & TREE-QMC-n2 vs. wASTRID & & 0.7 & & 43 & 3 & 4 & & 2e-11 & ***** & MC &  \\
FNR & TREE-QMC-n2 vs. wASTRID & & 0.75 & & 46 & 1 & 3 & & 7e-14 & ***** & MC &  \\

FPR & TREE-QMC-n2 vs. wASTRID & & 0.5 & & 31 & 10 & 9 & & 2e-04 & *** &  &  \\
FPR & TREE-QMC-n2 vs. wASTRID & & 0.55 & & 33 & 15 & 2 & & 9e-04 & ** &  &  \\
FPR & TREE-QMC-n2 vs. wASTRID & & 0.6 & & 39 & 9 & 2 & & 1e-08 & ***** & MC &  \\
FPR & TREE-QMC-n2 vs. wASTRID & & 0.65 & & 45 & 5 & 0 & & 4e-13 & ***** & MC &  \\
FPR & TREE-QMC-n2 vs. wASTRID & & 0.7 & & 49 & 1 & 0 & & 4e-15 & ***** & MC &  \\
FPR & TREE-QMC-n2 vs. wASTRID & & 0.75 & & 49 & 1 & 0 & & 4e-15 & ***** & MC &  \\


FNR & TREE-QMC-n2 vs. TREE-QMC-n2-shared & & 0.5 & & 34 & 7 & 9 & & 2e-06 & ***** & MC &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n2-shared & & 0.55 & & 39 & 5 & 6 & & 1e-09 & ***** & MC &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n2-shared & & 0.6 & & 37 & 7 & 6 & & 5e-08 & ***** & MC &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n2-shared & & 0.65 & & 45 & 5 & 0 & & 1e-10 & ***** & MC &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n2-shared & & 0.7 & & 44 & 1 & 5 & & 2e-13 & ***** & MC &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n2-shared & & 0.75 & & 41 & 6 & 3 & & 4e-09 & ***** & MC &  \\

FPR & TREE-QMC-n2 vs. TREE-QMC-n2-shared & & 0.5 & & 34 & 7 & 9 & & 3e-06 & ***** & MC &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n2-shared & & 0.55 & & 38 & 6 & 6 & & 2e-09 & ***** & MC &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n2-shared & & 0.6 & & 39 & 7 & 4 & & 1e-08 & ***** & MC &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n2-shared & & 0.65 & & 46 & 4 & 0 & & 1e-11 & ***** & MC &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n2-shared & & 0.7 & & 46 & 3 & 1 & & 5e-13 & ***** & MC &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n2-shared & & 0.75 & & 44 & 6 & 0 & & 2e-10 & ***** & MC &  \\


FNR & TREE-QMC-n2 vs. TREE-QMC-n1 & & 0.5 & & 16 & 12 & 22 & & 0.4 &  &  &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n1 & & 0.55 & & 16 & 12 & 22 & & 0.7 &  &  &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n1 & & 0.6 & & 11 & 13 & 26 & & 0.9 &  &  &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n1 & & 0.65 & & 21 & 12 & 17 & & 0.09 &  &  &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n1 & & 0.7 & & 21 & 15 & 14 & & 0.3 &  &  &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n1 & & 0.75 & & 14 & 22 & 14 & & 0.2 &  &  &  \\

FPR & TREE-QMC-n2 vs. TREE-QMC-n1 & & 0.5 & & 16 & 12 & 22 & & 0.5 &  &  &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n1 & & 0.55 & & 17 & 12 & 21 & & 0.8 &  &  &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n1 & & 0.6 & & 13 & 13 & 24 & & 1 &  &  &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n1 & & 0.65 & & 22 & 12 & 16 & & 0.05 &  &  &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n1 & & 0.7 & & 24 & 15 & 11 & & 0.3 &  &  &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n1 & & 0.75 & & 18 & 24 & 8 & & 0.3 &  &  &  \\


FNR & TREE-QMC-n2 vs. TREE-QMC-n0 & & 0.5 & & 33 & 8 & 9 & & 0.001 & ** &  &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n0 & & 0.55 & & 31 & 6 & 13 & & 3e-05 & **** & MC &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n0 & & 0.6 & & 30 & 12 & 8 & & 6e-04 & ** &  &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n0 & & 0.65 & & 38 & 8 & 4 & & 3e-07 & ***** & MC &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n0 & & 0.7 & & 31 & 11 & 8 & & 0.003 & ** &  &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n0 & & 0.75 & & 30 & 15 & 5 & & 0.01 & * &  &  \\

FPR & TREE-QMC-n2 vs. TREE-QMC-n0 & & 0.5 & & 33 & 9 & 8 & & 0.002 & ** &  &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n0 & & 0.55 & & 35 & 6 & 9 & & 1e-05 & **** & MC &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n0 & & 0.6 & & 32 & 13 & 5 & & 5e-04 & ** &  &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n0 & & 0.65 & & 40 & 7 & 3 & & 2e-07 & ***** & MC &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n0 & & 0.7 & & 32 & 13 & 5 & & 9e-04 & ** &  &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n0 & & 0.75 & & 36 & 13 & 1 & & 5e-04 & ** &  &  \\



