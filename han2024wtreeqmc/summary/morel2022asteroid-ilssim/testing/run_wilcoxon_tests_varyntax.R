#!/usr/bin/env Rscript

library(exactRankTests)
library(coin)

data <- read.csv("../csvs/data-varyntax-for-testing.csv")
data$NTAX <- as.factor(data$NTAX)
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
        for (ntax in ntaxs) {
        
        xdf <- data[data$NTAX == ntax, ]

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
                         as.character(ntax), "& &",
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
FNR & TREE-QMC-n2 vs. ASTEROID & & 25 & & 13 & 18 & 19 & & 0.3 &  &  &  \\
FNR & TREE-QMC-n2 vs. ASTEROID & & 75 & & 14 & 27 & 9 & & 0.008 & * &  & ( ASTEROID better) \\
FNR & TREE-QMC-n2 vs. ASTEROID & & 50 & & 12 & 27 & 11 & & 0.01 & * &  & ( ASTEROID better) \\
FNR & TREE-QMC-n2 vs. ASTEROID & & 100 & & 13 & 28 & 9 & & 0.1 &  &  &  \\
FNR & TREE-QMC-n2 vs. ASTEROID & & 125 & & 13 & 29 & 8 & & 0.003 & ** &  & ( ASTEROID better) \\
FNR & TREE-QMC-n2 vs. ASTEROID & & 150 & & 15 & 32 & 3 & & 0.002 & ** &  & ( ASTEROID better) \\

FPR & TREE-QMC-n2 vs. ASTEROID & & 25 & & 20 & 16 & 14 & & 0.7 &  &  &  \\
FPR & TREE-QMC-n2 vs. ASTEROID & & 75 & & 23 & 25 & 2 & & 0.8 &  &  &  \\
FPR & TREE-QMC-n2 vs. ASTEROID & & 50 & & 24 & 23 & 3 & & 0.7 &  &  &  \\
FPR & TREE-QMC-n2 vs. ASTEROID & & 100 & & 29 & 20 & 1 & & 0.03 & * &  &  \\
FPR & TREE-QMC-n2 vs. ASTEROID & & 125 & & 30 & 19 & 1 & & 0.01 & * &  &  \\
FPR & TREE-QMC-n2 vs. ASTEROID & & 150 & & 31 & 19 & 0 & & 0.1 &  &  &  \\


FNR & TREE-QMC-n2 vs. ASTER & & 25 & & 29 & 9 & 12 & & 8e-04 & ** &  &  \\
FNR & TREE-QMC-n2 vs. ASTER & & 75 & & 26 & 18 & 6 & & 0.1 &  &  &  \\
FNR & TREE-QMC-n2 vs. ASTER & & 50 & & 24 & 19 & 7 & & 0.3 &  &  &  \\
FNR & TREE-QMC-n2 vs. ASTER & & 100 & & 30 & 15 & 5 & & 0.01 & * &  &  \\
FNR & TREE-QMC-n2 vs. ASTER & & 125 & & 31 & 16 & 3 & & 0.003 & ** &  &  \\
FNR & TREE-QMC-n2 vs. ASTER & & 150 & & 38 & 6 & 6 & & 3e-06 & ***** & MC &  \\

FPR & TREE-QMC-n2 vs. ASTER & & 25 & & 36 & 7 & 7 & & 8e-06 & **** & MC &  \\
FPR & TREE-QMC-n2 vs. ASTER & & 75 & & 33 & 16 & 1 & & 3e-04 & *** &  &  \\
FPR & TREE-QMC-n2 vs. ASTER & & 50 & & 28 & 18 & 4 & & 0.009 & * &  &  \\
FPR & TREE-QMC-n2 vs. ASTER & & 100 & & 39 & 10 & 1 & & 4e-08 & ***** & MC &  \\
FPR & TREE-QMC-n2 vs. ASTER & & 125 & & 41 & 9 & 0 & & 2e-08 & ***** & MC &  \\
FPR & TREE-QMC-n2 vs. ASTER & & 150 & & 46 & 4 & 0 & & 4e-13 & ***** & MC &  \\


FNR & TREE-QMC-n2 vs. wASTRID & & 25 & & 25 & 14 & 11 & & 0.04 & * &  &  \\
FNR & TREE-QMC-n2 vs. wASTRID & & 75 & & 36 & 9 & 5 & & 5e-07 & ***** & MC &  \\
FNR & TREE-QMC-n2 vs. wASTRID & & 50 & & 35 & 10 & 5 & & 1e-06 & ***** & MC &  \\
FNR & TREE-QMC-n2 vs. wASTRID & & 100 & & 40 & 8 & 2 & & 6e-09 & ***** & MC &  \\
FNR & TREE-QMC-n2 vs. wASTRID & & 125 & & 45 & 4 & 1 & & 3e-11 & ***** & MC &  \\
FNR & TREE-QMC-n2 vs. wASTRID & & 150 & & 48 & 2 & 0 & & 5e-13 & ***** & MC &  \\

FPR & TREE-QMC-n2 vs. wASTRID & & 25 & & 31 & 12 & 7 & & 0.002 & ** &  &  \\
FPR & TREE-QMC-n2 vs. wASTRID & & 75 & & 42 & 7 & 1 & & 9e-10 & ***** & MC &  \\
FPR & TREE-QMC-n2 vs. wASTRID & & 50 & & 39 & 9 & 2 & & 1e-08 & ***** & MC &  \\
FPR & TREE-QMC-n2 vs. wASTRID & & 100 & & 43 & 6 & 1 & & 1e-11 & ***** & MC &  \\
FPR & TREE-QMC-n2 vs. wASTRID & & 125 & & 47 & 3 & 0 & & 2e-13 & ***** & MC &  \\
FPR & TREE-QMC-n2 vs. wASTRID & & 150 & & 49 & 1 & 0 & & 9e-15 & ***** & MC &  \\


FNR & TREE-QMC-n2 vs. TREE-QMC-n2-shared & & 25 & & 28 & 9 & 13 & & 1e-04 & *** & MC &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n2-shared & & 75 & & 46 & 2 & 2 & & 1e-12 & ***** & MC &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n2-shared & & 50 & & 37 & 7 & 6 & & 5e-08 & ***** & MC &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n2-shared & & 100 & & 50 & 0 & 0 & & 2e-15 & ***** & MC &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n2-shared & & 125 & & 50 & 0 & 0 & & 2e-15 & ***** & MC &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n2-shared & & 150 & & 50 & 0 & 0 & & 2e-15 & ***** & MC &  \\

FPR & TREE-QMC-n2 vs. TREE-QMC-n2-shared & & 25 & & 28 & 11 & 11 & & 2e-04 & *** &  &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n2-shared & & 75 & & 46 & 3 & 1 & & 7e-13 & ***** & MC &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n2-shared & & 50 & & 39 & 7 & 4 & & 1e-08 & ***** & MC &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n2-shared & & 100 & & 50 & 0 & 0 & & 2e-15 & ***** & MC &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n2-shared & & 125 & & 50 & 0 & 0 & & 2e-15 & ***** & MC &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n2-shared & & 150 & & 50 & 0 & 0 & & 2e-15 & ***** & MC &  \\


FNR & TREE-QMC-n2 vs. TREE-QMC-n1 & & 25 & & 7 & 8 & 35 & & 1 &  &  &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n1 & & 75 & & 22 & 15 & 13 & & 0.2 &  &  &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n1 & & 50 & & 11 & 13 & 26 & & 0.9 &  &  &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n1 & & 100 & & 21 & 17 & 12 & & 0.2 &  &  &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n1 & & 125 & & 26 & 16 & 8 & & 0.04 & * &  &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n1 & & 150 & & 27 & 16 & 7 & & 0.03 & * &  &  \\

FPR & TREE-QMC-n2 vs. TREE-QMC-n1 & & 25 & & 7 & 8 & 35 & & 1 &  &  &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n1 & & 75 & & 22 & 16 & 12 & & 0.3 &  &  &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n1 & & 50 & & 13 & 13 & 24 & & 1 &  &  &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n1 & & 100 & & 21 & 17 & 12 & & 0.2 &  &  &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n1 & & 125 & & 27 & 16 & 7 & & 0.2 &  &  &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n1 & & 150 & & 28 & 16 & 6 & & 0.03 & * &  &  \\


FNR & TREE-QMC-n2 vs. TREE-QMC-n0 & & 25 & & 29 & 8 & 13 & & 6e-05 & *** & MC &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n0 & & 75 & & 30 & 11 & 9 & & 3e-04 & *** &  &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n0 & & 50 & & 30 & 12 & 8 & & 6e-04 & ** &  &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n0 & & 100 & & 39 & 7 & 4 & & 7e-09 & ***** & MC &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n0 & & 125 & & 39 & 6 & 5 & & 1e-07 & ***** & MC &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n0 & & 150 & & 46 & 3 & 1 & & 7e-13 & ***** & MC &  \\

FPR & TREE-QMC-n2 vs. TREE-QMC-n0 & & 25 & & 30 & 10 & 10 & & 3e-04 & *** &  &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n0 & & 75 & & 30 & 13 & 7 & & 4e-04 & *** &  &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n0 & & 50 & & 32 & 13 & 5 & & 5e-04 & ** &  &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n0 & & 100 & & 40 & 8 & 2 & & 3e-09 & ***** & MC &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n0 & & 125 & & 42 & 5 & 3 & & 5e-08 & ***** & MC &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n0 & & 150 & & 47 & 3 & 0 & & 8e-14 & ***** & MC &  \\


[1] "done"

