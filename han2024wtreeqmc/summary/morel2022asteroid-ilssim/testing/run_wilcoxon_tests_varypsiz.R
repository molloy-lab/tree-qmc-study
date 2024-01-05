#!/usr/bin/env Rscript

library(exactRankTests)
library(coin)

data <- read.csv("../csvs/data-varypsiz-for-testing.csv")
data$PSIZ <- as.factor(data$PSIZ)
data$REPL <- as.factor(data$REPL)

print(data)

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

        writeme <- paste("TQMC-wh-n2 vs.", mthd, "&",
                         emet, "&",
                         as.character(psiz), "&",
                         as.character(ntreeqmc), "/",
                         as.character(nother), "/",
                         as.character(ntie), "&",
                         format(xwsr, scientific = TRUE), "&", 
                         format(ywsr, scientific = TRUE), "&",
                         "&", stars, "&", mc, note)
        print(writeme)

        }
        print("")
    }
    print("")
}

print("done")
exit()

[1] "alpha = 0.05"
[1] "bonferroni alpha = 0.00016025641025641 = 0.05 / 312"
[1] "TQMC-wh-n2 vs. ASTEROID & FNR & 10 & 15 / 24 / 11 & 0.02719862 & 0.04426908 & & * &  ( ASTEROID better)"
[1] "TQMC-wh-n2 vs. ASTEROID & FNR & 50000000 & 12 / 27 / 11 & 0.01183831 & 0.01055658 & & * &  ( ASTEROID better)"
[1] "TQMC-wh-n2 vs. ASTEROID & FNR & 100000000 & 20 / 20 / 10 & 0.2536306 & 0.4111414 & &  &  "
[1] "TQMC-wh-n2 vs. ASTEROID & FNR & 500000000 & 26 / 20 / 4 & 0.2569219 & 0.2637694 & &  &  "
[1] "TQMC-wh-n2 vs. ASTEROID & FNR & 1000000000 & 16 / 25 / 9 & 0.3728778 & 0.2902527 & &  &  "
[1] ""
[1] "TQMC-wh-n2 vs. ASTEROID & FPR & 10 & 20 / 23 / 7 & 0.9785477 & 0.9027547 & &  &  "
[1] "TQMC-wh-n2 vs. ASTEROID & FPR & 50000000 & 24 / 23 / 3 & 0.6690783 & 0.7071011 & &  &  "
[1] "TQMC-wh-n2 vs. ASTEROID & FPR & 100000000 & 30 / 19 / 1 & 0.1751631 & 0.1710362 & &  &  "
[1] "TQMC-wh-n2 vs. ASTEROID & FPR & 500000000 & 30 / 18 / 2 & 0.03671343 & 0.03734623 & & * &  "
[1] "TQMC-wh-n2 vs. ASTEROID & FPR & 1000000000 & 22 / 24 / 4 & 0.8307513 & 0.8201525 & &  &  "
[1] ""
[1] ""
[1] "TQMC-wh-n2 vs. ASTER & FNR & 10 & 30 / 13 / 7 & 0.005180817 & 0.004836959 & & * &  "
[1] "TQMC-wh-n2 vs. ASTER & FNR & 50000000 & 24 / 19 / 7 & 0.3086718 & 0.3253 & &  &  "
[1] "TQMC-wh-n2 vs. ASTER & FNR & 100000000 & 26 / 16 / 8 & 0.01304002 & 0.02138751 & & * &  "
[1] "TQMC-wh-n2 vs. ASTER & FNR & 500000000 & 32 / 12 / 6 & 0.0004457919 & 0.0004816404 & & *** &  "
[1] "TQMC-wh-n2 vs. ASTER & FNR & 1000000000 & 32 / 16 / 2 & 0.005113407 & 0.00527816 & & * &  "
[1] ""
[1] "TQMC-wh-n2 vs. ASTER & FPR & 10 & 36 / 12 / 2 & 2.760526e-05 & 2.928365e-05 & & **** & MC "
[1] "TQMC-wh-n2 vs. ASTER & FPR & 50000000 & 28 / 18 / 4 & 0.006314081 & 0.009292476 & & * &  "
[1] "TQMC-wh-n2 vs. ASTER & FPR & 100000000 & 36 / 14 / 0 & 5.167958e-05 & 5.167958e-05 & & *** & MC "
[1] "TQMC-wh-n2 vs. ASTER & FPR & 500000000 & 40 / 10 / 0 & 1.208874e-05 & 1.208874e-05 & & **** & MC "
[1] "TQMC-wh-n2 vs. ASTER & FPR & 1000000000 & 32 / 16 / 2 & 0.0001553798 & 0.0002091479 & & *** &  "
[1] ""
[1] ""
[1] "TQMC-wh-n2 vs. wASTRID & FNR & 10 & 40 / 5 / 5 & 2.901572e-09 & 2.550223e-09 & & ***** & MC "
[1] "TQMC-wh-n2 vs. wASTRID & FNR & 50000000 & 35 / 10 / 5 & 6.915856e-07 & 1.321502e-06 & & ***** & MC "
[1] "TQMC-wh-n2 vs. wASTRID & FNR & 100000000 & 38 / 8 / 4 & 1.123395e-06 & 9.548156e-07 & & ***** & MC "
[1] "TQMC-wh-n2 vs. wASTRID & FNR & 500000000 & 34 / 11 / 5 & 4.100076e-05 & 4.726777e-05 & & **** & MC "
[1] "TQMC-wh-n2 vs. wASTRID & FNR & 1000000000 & 30 / 15 / 5 & 0.053088 & 0.04340889 & &  &  "
[1] ""
[1] "TQMC-wh-n2 vs. wASTRID & FPR & 10 & 44 / 4 / 2 & 6.103562e-12 & 6.743051e-12 & & ***** & MC "
[1] "TQMC-wh-n2 vs. wASTRID & FPR & 50000000 & 39 / 9 / 2 & 7.50876e-09 & 1.083473e-08 & & ***** & MC "
[1] "TQMC-wh-n2 vs. wASTRID & FPR & 100000000 & 40 / 8 / 2 & 2.278426e-08 & 2.507222e-08 & & ***** & MC "
[1] "TQMC-wh-n2 vs. wASTRID & FPR & 500000000 & 38 / 10 / 2 & 5.035047e-06 & 4.841724e-06 & & **** & MC "
[1] "TQMC-wh-n2 vs. wASTRID & FPR & 1000000000 & 34 / 15 / 1 & 0.006106057 & 0.005891067 & & * &  "
[1] ""
[1] ""
[1] "TQMC-wh-n2 vs. TREE-QMC-n2-shared & FNR & 10 & 41 / 5 / 4 & 1.394767e-09 & 1.228443e-09 & & ***** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n2-shared & FNR & 50000000 & 37 / 7 / 6 & 2.857598e-08 & 4.626281e-08 & & ***** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n2-shared & FNR & 100000000 & 45 / 1 / 4 & 1.98952e-13 & 1.421085e-13 & & ***** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n2-shared & FNR & 500000000 & 44 / 3 / 3 & 1.545004e-10 & 8.452616e-11 & & ***** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n2-shared & FNR & 1000000000 & 46 / 4 / 0 & 4.312106e-11 & 4.312106e-11 & & ***** & MC "
[1] ""
[1] "TQMC-wh-n2 vs. TREE-QMC-n2-shared & FPR & 10 & 42 / 5 / 3 & 7.025704e-10 & 6.254055e-10 & & ***** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n2-shared & FPR & 50000000 & 39 / 7 / 4 & 7.394078e-09 & 1.012964e-08 & & ***** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n2-shared & FPR & 100000000 & 46 / 1 / 3 & 4.263256e-14 & 4.263256e-14 & & ***** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n2-shared & FPR & 500000000 & 45 / 4 / 1 & 1.978826e-10 & 1.658123e-10 & & ***** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n2-shared & FPR & 1000000000 & 46 / 4 / 0 & 2.917311e-11 & 2.917311e-11 & & ***** & MC "
[1] ""
[1] ""
[1] "TQMC-wh-n2 vs. TREE-QMC-n1 & FNR & 10 & 17 / 12 / 21 & 0.5884478 & 0.4380318 & &  &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n1 & FNR & 50000000 & 11 / 13 / 26 & 0.903972 & 0.8266549 & &  &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n1 & FNR & 100000000 & 18 / 16 / 16 & 0.4686499 & 0.5667434 & &  &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n1 & FNR & 500000000 & 19 / 17 / 14 & 0.3189716 & 0.4503818 & &  &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n1 & FNR & 1000000000 & 14 / 15 / 21 & 0.7646109 & 0.8072375 & &  &  "
[1] ""
[1] "TQMC-wh-n2 vs. TREE-QMC-n1 & FPR & 10 & 17 / 13 / 20 & 0.499097 & 0.4678225 & &  &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n1 & FPR & 50000000 & 13 / 13 / 24 & 0.8950852 & 0.9571749 & &  &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n1 & FPR & 100000000 & 18 / 17 / 15 & 0.6532241 & 0.7321901 & &  &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n1 & FPR & 500000000 & 19 / 17 / 14 & 0.4479936 & 0.5430702 & &  &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n1 & FPR & 1000000000 & 17 / 18 / 15 & 0.8808867 & 0.8707998 & &  &  "
[1] ""
[1] ""
[1] "TQMC-wh-n2 vs. TREE-QMC-n0 & FNR & 10 & 33 / 10 / 7 & 3.366892e-05 & 4.080925e-05 & & **** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n0 & FNR & 50000000 & 30 / 12 / 8 & 0.0003918257 & 0.0005882718 & & ** &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n0 & FNR & 100000000 & 33 / 11 / 6 & 1.513519e-05 & 2.519368e-05 & & **** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n0 & FNR & 500000000 & 33 / 8 / 9 & 0.0001200018 & 6.371319e-05 & & *** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n0 & FNR & 1000000000 & 35 / 11 / 4 & 2.100755e-06 & 3.400438e-06 & & ***** & MC "
[1] ""
[1] "TQMC-wh-n2 vs. TREE-QMC-n0 & FPR & 10 & 35 / 9 / 6 & 2.701806e-05 & 2.151234e-05 & & **** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n0 & FPR & 50000000 & 32 / 13 / 5 & 0.0004205377 & 0.0005079325 & & ** &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n0 & FPR & 100000000 & 35 / 11 / 4 & 1.175277e-05 & 1.444829e-05 & & **** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n0 & FPR & 500000000 & 34 / 11 / 5 & 0.0001268964 & 0.0001214938 & & *** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n0 & FPR & 1000000000 & 36 / 14 / 0 & 2.260878e-06 & 2.260878e-06 & & ***** & MC "
[1] ""
[1] ""
[1] "done"
