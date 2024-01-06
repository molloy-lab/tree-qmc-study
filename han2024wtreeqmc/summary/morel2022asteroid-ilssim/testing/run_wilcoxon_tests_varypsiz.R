#!/usr/bin/env Rscript

library(exactRankTests)
library(coin)

data <- read.csv("../csvs/data-varypsiz-for-testing.csv")
data$PSIZ <- as.factor(data$PSIZ)
data$REPL <- as.factor(data$REPL)

mthds <- c("ASTEROID", "ASTER", "wASTRID",
           "TREE-QMC-n2-shared", "TREE-QMC-n1", "TREE-QMC-n0")
emets <- c("FNR", "FPR")

psizs <- c("10", "50000000", "100000000", "500000000", "1000000000")  # BASE CONDITION: 50000000
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
        for (psiz in psizs) {
        
        xdf <- data[data$PSIZ == psiz, ]

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
                         as.character(psiz), "& &",
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
FNR & TREE-QMC-n2 vs. ASTEROID & & 10 & & 15 & 24 & 11 & & 0.04 & * &  & ( ASTEROID better) \\
FNR & TREE-QMC-n2 vs. ASTEROID & & 50000000 & & 12 & 27 & 11 & & 0.01 & * &  & ( ASTEROID better) \\
FNR & TREE-QMC-n2 vs. ASTEROID & & 100000000 & & 20 & 20 & 10 & & 0.4 &  &  &  \\
FNR & TREE-QMC-n2 vs. ASTEROID & & 500000000 & & 26 & 20 & 4 & & 0.3 &  &  &  \\
FNR & TREE-QMC-n2 vs. ASTEROID & & 1000000000 & & 16 & 25 & 9 & & 0.4 &  &  &  \\

FPR & TREE-QMC-n2 vs. ASTEROID & & 10 & & 20 & 23 & 7 & & 1 &  &  &  \\
FPR & TREE-QMC-n2 vs. ASTEROID & & 50000000 & & 24 & 23 & 3 & & 0.7 &  &  &  \\
FPR & TREE-QMC-n2 vs. ASTEROID & & 100000000 & & 30 & 19 & 1 & & 0.2 &  &  &  \\
FPR & TREE-QMC-n2 vs. ASTEROID & & 500000000 & & 30 & 18 & 2 & & 0.04 & * &  &  \\
FPR & TREE-QMC-n2 vs. ASTEROID & & 1000000000 & & 22 & 24 & 4 & & 0.8 &  &  &  \\


FNR & TREE-QMC-n2 vs. ASTER & & 10 & & 30 & 13 & 7 & & 0.005 & * &  &  \\
FNR & TREE-QMC-n2 vs. ASTER & & 50000000 & & 24 & 19 & 7 & & 0.3 &  &  &  \\
FNR & TREE-QMC-n2 vs. ASTER & & 100000000 & & 26 & 16 & 8 & & 0.02 & * &  &  \\
FNR & TREE-QMC-n2 vs. ASTER & & 500000000 & & 32 & 12 & 6 & & 5e-04 & *** &  &  \\
FNR & TREE-QMC-n2 vs. ASTER & & 1000000000 & & 32 & 16 & 2 & & 0.005 & * &  &  \\

FPR & TREE-QMC-n2 vs. ASTER & & 10 & & 36 & 12 & 2 & & 3e-05 & **** & MC &  \\
FPR & TREE-QMC-n2 vs. ASTER & & 50000000 & & 28 & 18 & 4 & & 0.009 & * &  &  \\
FPR & TREE-QMC-n2 vs. ASTER & & 100000000 & & 36 & 14 & 0 & & 5e-05 & *** & MC &  \\
FPR & TREE-QMC-n2 vs. ASTER & & 500000000 & & 40 & 10 & 0 & & 1e-05 & **** & MC &  \\
FPR & TREE-QMC-n2 vs. ASTER & & 1000000000 & & 32 & 16 & 2 & & 2e-04 & *** &  &  \\


FNR & TREE-QMC-n2 vs. wASTRID & & 10 & & 40 & 5 & 5 & & 3e-09 & ***** & MC &  \\
FNR & TREE-QMC-n2 vs. wASTRID & & 50000000 & & 35 & 10 & 5 & & 1e-06 & ***** & MC &  \\
FNR & TREE-QMC-n2 vs. wASTRID & & 100000000 & & 38 & 8 & 4 & & 1e-06 & ***** & MC &  \\
FNR & TREE-QMC-n2 vs. wASTRID & & 500000000 & & 34 & 11 & 5 & & 5e-05 & **** & MC &  \\
FNR & TREE-QMC-n2 vs. wASTRID & & 1000000000 & & 30 & 15 & 5 & & 0.05 &  &  &  \\

FPR & TREE-QMC-n2 vs. wASTRID & & 10 & & 44 & 4 & 2 & & 7e-12 & ***** & MC &  \\
FPR & TREE-QMC-n2 vs. wASTRID & & 50000000 & & 39 & 9 & 2 & & 1e-08 & ***** & MC &  \\
FPR & TREE-QMC-n2 vs. wASTRID & & 100000000 & & 40 & 8 & 2 & & 3e-08 & ***** & MC &  \\
FPR & TREE-QMC-n2 vs. wASTRID & & 500000000 & & 38 & 10 & 2 & & 5e-06 & **** & MC &  \\
FPR & TREE-QMC-n2 vs. wASTRID & & 1000000000 & & 34 & 15 & 1 & & 0.006 & * &  &  \\


FNR & TREE-QMC-n2 vs. TREE-QMC-n2-shared & & 10 & & 41 & 5 & 4 & & 1e-09 & ***** & MC &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n2-shared & & 50000000 & & 37 & 7 & 6 & & 5e-08 & ***** & MC &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n2-shared & & 100000000 & & 45 & 1 & 4 & & 2e-13 & ***** & MC &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n2-shared & & 500000000 & & 44 & 3 & 3 & & 2e-10 & ***** & MC &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n2-shared & & 1000000000 & & 46 & 4 & 0 & & 4e-11 & ***** & MC &  \\

FPR & TREE-QMC-n2 vs. TREE-QMC-n2-shared & & 10 & & 42 & 5 & 3 & & 7e-10 & ***** & MC &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n2-shared & & 50000000 & & 39 & 7 & 4 & & 1e-08 & ***** & MC &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n2-shared & & 100000000 & & 46 & 1 & 3 & & 4e-14 & ***** & MC &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n2-shared & & 500000000 & & 45 & 4 & 1 & & 2e-10 & ***** & MC &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n2-shared & & 1000000000 & & 46 & 4 & 0 & & 3e-11 & ***** & MC &  \\


FNR & TREE-QMC-n2 vs. TREE-QMC-n1 & & 10 & & 17 & 12 & 21 & & 0.6 &  &  &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n1 & & 50000000 & & 11 & 13 & 26 & & 0.9 &  &  &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n1 & & 100000000 & & 18 & 16 & 16 & & 0.6 &  &  &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n1 & & 500000000 & & 19 & 17 & 14 & & 0.5 &  &  &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n1 & & 1000000000 & & 14 & 15 & 21 & & 0.8 &  &  &  \\

FPR & TREE-QMC-n2 vs. TREE-QMC-n1 & & 10 & & 17 & 13 & 20 & & 0.5 &  &  &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n1 & & 50000000 & & 13 & 13 & 24 & & 1 &  &  &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n1 & & 100000000 & & 18 & 17 & 15 & & 0.7 &  &  &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n1 & & 500000000 & & 19 & 17 & 14 & & 0.5 &  &  &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n1 & & 1000000000 & & 17 & 18 & 15 & & 0.9 &  &  &  \\


FNR & TREE-QMC-n2 vs. TREE-QMC-n0 & & 10 & & 33 & 10 & 7 & & 4e-05 & **** & MC &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n0 & & 50000000 & & 30 & 12 & 8 & & 6e-04 & ** &  &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n0 & & 100000000 & & 33 & 11 & 6 & & 3e-05 & **** & MC &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n0 & & 500000000 & & 33 & 8 & 9 & & 1e-04 & *** & MC &  \\
FNR & TREE-QMC-n2 vs. TREE-QMC-n0 & & 1000000000 & & 35 & 11 & 4 & & 3e-06 & ***** & MC &  \\

FPR & TREE-QMC-n2 vs. TREE-QMC-n0 & & 10 & & 35 & 9 & 6 & & 3e-05 & **** & MC &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n0 & & 50000000 & & 32 & 13 & 5 & & 5e-04 & ** &  &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n0 & & 100000000 & & 35 & 11 & 4 & & 1e-05 & **** & MC &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n0 & & 500000000 & & 34 & 11 & 5 & & 1e-04 & *** & MC &  \\
FPR & TREE-QMC-n2 vs. TREE-QMC-n0 & & 1000000000 & & 36 & 14 & 0 & & 2e-06 & ***** & MC &  \\


[1] "done"

