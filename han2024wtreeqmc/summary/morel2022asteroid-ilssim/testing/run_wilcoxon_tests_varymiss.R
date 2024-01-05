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

        writeme <- paste("TQMC-wh-n2 vs.", mthd, "&",
                         emet, "&",
                         as.character(miss), "&",
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
[1] "TQMC-wh-n2 vs. ASTEROID & FNR & 0.5 & 14 / 23 / 13 & 0.2575892 & 0.1932711 & &  &  "
[1] "TQMC-wh-n2 vs. ASTEROID & FNR & 0.55 & 15 / 17 / 18 & 0.7837584 & 0.7474729 & &  &  "
[1] "TQMC-wh-n2 vs. ASTEROID & FNR & 0.6 & 12 / 27 / 11 & 0.01183831 & 0.01055658 & & * &  ( ASTEROID better)"
[1] "TQMC-wh-n2 vs. ASTEROID & FNR & 0.65 & 15 / 24 / 11 & 0.07283367 & 0.08412945 & &  &  "
[1] "TQMC-wh-n2 vs. ASTEROID & FNR & 0.7 & 10 / 32 / 8 & 1.12646e-06 & 4.599623e-06 & & ***** & MC ( ASTEROID better)"
[1] "TQMC-wh-n2 vs. ASTEROID & FNR & 0.75 & 7 / 40 / 3 & 3.851497e-09 & 4.828934e-09 & & ***** & MC ( ASTEROID better)"
[1] ""
[1] "TQMC-wh-n2 vs. ASTEROID & FPR & 0.5 & 17 / 23 / 10 & 0.2827352 & 0.2870921 & &  &  "
[1] "TQMC-wh-n2 vs. ASTEROID & FPR & 0.55 & 22 / 16 / 12 & 0.271665 & 0.2773415 & &  &  "
[1] "TQMC-wh-n2 vs. ASTEROID & FPR & 0.6 & 24 / 23 / 3 & 0.6690783 & 0.7071011 & &  &  "
[1] "TQMC-wh-n2 vs. ASTEROID & FPR & 0.65 & 39 / 11 / 0 & 1.031201e-05 & 1.031201e-05 & & **** & MC "
[1] "TQMC-wh-n2 vs. ASTEROID & FPR & 0.7 & 33 / 17 / 0 & 0.07237748 & 0.07237748 & &  &  "
[1] "TQMC-wh-n2 vs. ASTEROID & FPR & 0.75 & 28 / 22 / 0 & 0.153598 & 0.153598 & &  &  "
[1] ""
[1] ""
[1] "TQMC-wh-n2 vs. ASTER & FNR & 0.5 & 27 / 11 / 12 & 0.01632847 & 0.01067048 & & * &  "
[1] "TQMC-wh-n2 vs. ASTER & FNR & 0.55 & 31 / 9 / 10 & 0.0004721495 & 0.0003047211 & & *** &  "
[1] "TQMC-wh-n2 vs. ASTER & FNR & 0.6 & 24 / 19 / 7 & 0.3086718 & 0.3253 & &  &  "
[1] "TQMC-wh-n2 vs. ASTER & FNR & 0.65 & 29 / 7 / 14 & 0.0001563749 & 9.944837e-05 & & *** & MC "
[1] "TQMC-wh-n2 vs. ASTER & FNR & 0.7 & 22 / 20 / 8 & 0.376453 & 0.4463233 & &  &  "
[1] "TQMC-wh-n2 vs. ASTER & FNR & 0.75 & 28 / 15 / 7 & 0.00788653 & 0.0103285 & & * &  "
[1] ""
[1] "TQMC-wh-n2 vs. ASTER & FPR & 0.5 & 28 / 11 / 11 & 0.009460015 & 0.006461527 & & * &  "
[1] "TQMC-wh-n2 vs. ASTER & FPR & 0.55 & 36 / 9 / 5 & 2.76519e-05 & 2.054107e-05 & & **** & MC "
[1] "TQMC-wh-n2 vs. ASTER & FPR & 0.6 & 28 / 18 / 4 & 0.006314081 & 0.009292476 & & * &  "
[1] "TQMC-wh-n2 vs. ASTER & FPR & 0.65 & 45 / 4 / 1 & 1.257661e-12 & 1.41398e-12 & & ***** & MC "
[1] "TQMC-wh-n2 vs. ASTER & FPR & 0.7 & 43 / 7 / 0 & 7.546861e-09 & 7.546861e-09 & & ***** & MC "
[1] "TQMC-wh-n2 vs. ASTER & FPR & 0.75 & 45 / 5 / 0 & 3.487575e-10 & 3.487575e-10 & & ***** & MC "
[1] ""
[1] ""
[1] "TQMC-wh-n2 vs. wASTRID & FNR & 0.5 & 29 / 10 / 11 & 0.000280377 & 0.0003921684 & & *** &  "
[1] "TQMC-wh-n2 vs. wASTRID & FNR & 0.55 & 32 / 15 / 3 & 0.003005621 & 0.003141185 & & ** &  "
[1] "TQMC-wh-n2 vs. wASTRID & FNR & 0.6 & 35 / 10 / 5 & 6.915856e-07 & 1.321502e-06 & & ***** & MC "
[1] "TQMC-wh-n2 vs. wASTRID & FNR & 0.65 & 35 / 9 / 6 & 6.554956e-09 & 4.160881e-08 & & ***** & MC "
[1] "TQMC-wh-n2 vs. wASTRID & FNR & 0.7 & 43 / 3 / 4 & 2.179945e-11 & 1.651301e-11 & & ***** & MC "
[1] "TQMC-wh-n2 vs. wASTRID & FNR & 0.75 & 46 / 1 / 3 & 7.105427e-14 & 7.105427e-14 & & ***** & MC "
[1] ""
[1] "TQMC-wh-n2 vs. wASTRID & FPR & 0.5 & 31 / 10 / 9 & 0.0002132821 & 0.0002221788 & & *** &  "
[1] "TQMC-wh-n2 vs. wASTRID & FPR & 0.55 & 33 / 15 / 2 & 0.0007973675 & 0.0008695461 & & ** &  "
[1] "TQMC-wh-n2 vs. wASTRID & FPR & 0.6 & 39 / 9 / 2 & 7.50876e-09 & 1.083473e-08 & & ***** & MC "
[1] "TQMC-wh-n2 vs. wASTRID & FPR & 0.65 & 45 / 5 / 0 & 4.423129e-13 & 4.423129e-13 & & ***** & MC "
[1] "TQMC-wh-n2 vs. wASTRID & FPR & 0.7 & 49 / 1 / 0 & 3.552714e-15 & 3.552714e-15 & & ***** & MC "
[1] "TQMC-wh-n2 vs. wASTRID & FPR & 0.75 & 49 / 1 / 0 & 3.552714e-15 & 3.552714e-15 & & ***** & MC "
[1] ""
[1] ""
[1] "TQMC-wh-n2 vs. TREE-QMC-n2-shared & FNR & 0.5 & 34 / 7 / 9 & 1.923996e-06 & 1.920574e-06 & & ***** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n2-shared & FNR & 0.55 & 39 / 5 / 6 & 6.035634e-10 & 9.630412e-10 & & ***** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n2-shared & FNR & 0.6 & 37 / 7 / 6 & 2.857598e-08 & 4.626281e-08 & & ***** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n2-shared & FNR & 0.65 & 45 / 5 / 0 & 1.02883e-10 & 1.02883e-10 & & ***** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n2-shared & FNR & 0.7 & 44 / 1 / 5 & 1.705303e-13 & 1.705303e-13 & & ***** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n2-shared & FNR & 0.75 & 41 / 6 / 3 & 4.264336e-09 & 4.040558e-09 & & ***** & MC "
[1] ""
[1] "TQMC-wh-n2 vs. TREE-QMC-n2-shared & FPR & 0.5 & 34 / 7 / 9 & 2.992451e-06 & 2.614127e-06 & & ***** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n2-shared & FPR & 0.55 & 38 / 6 / 6 & 6.190248e-10 & 1.789545e-09 & & ***** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n2-shared & FPR & 0.6 & 39 / 7 / 4 & 7.394078e-09 & 1.012964e-08 & & ***** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n2-shared & FPR & 0.65 & 46 / 4 / 0 & 1.052669e-11 & 1.052669e-11 & & ***** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n2-shared & FPR & 0.7 & 46 / 3 / 1 & 4.867218e-13 & 4.938272e-13 & & ***** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n2-shared & FPR & 0.75 & 44 / 6 / 0 & 1.615206e-10 & 1.615206e-10 & & ***** & MC "
[1] ""
[1] ""
[1] "TQMC-wh-n2 vs. TREE-QMC-n1 & FNR & 0.5 & 16 / 12 / 22 & 0.3689319 & 0.4000838 & &  &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n1 & FNR & 0.55 & 16 / 12 / 22 & 0.6977096 & 0.5392074 & &  &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n1 & FNR & 0.6 & 11 / 13 / 26 & 0.903972 & 0.8266549 & &  &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n1 & FNR & 0.65 & 21 / 12 / 17 & 0.08814617 & 0.09026367 & &  &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n1 & FNR & 0.7 & 21 / 15 / 14 & 0.2720685 & 0.2756909 & &  &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n1 & FNR & 0.75 & 14 / 22 / 14 & 0.2081196 & 0.1840792 & &  &  "
[1] ""
[1] "TQMC-wh-n2 vs. TREE-QMC-n1 & FPR & 0.5 & 16 / 12 / 22 & 0.5024127 & 0.4596193 & &  &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n1 & FPR & 0.55 & 17 / 12 / 21 & 0.8435464 & 0.5384593 & &  &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n1 & FPR & 0.6 & 13 / 13 / 24 & 0.8950852 & 0.9571749 & &  &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n1 & FPR & 0.65 & 22 / 12 / 16 & 0.04212491 & 0.05078951 & &  &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n1 & FPR & 0.7 & 24 / 15 / 11 & 0.2939408 & 0.2260485 & &  &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n1 & FPR & 0.75 & 18 / 24 / 8 & 0.2541658 & 0.2652408 & &  &  "
[1] ""
[1] ""
[1] "TQMC-wh-n2 vs. TREE-QMC-n0 & FNR & 0.5 & 33 / 8 / 9 & 0.001121154 & 0.0003697249 & & ** &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n0 & FNR & 0.55 & 31 / 6 / 13 & 2.513744e-05 & 1.295595e-05 & & **** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n0 & FNR & 0.6 & 30 / 12 / 8 & 0.0003918257 & 0.0005882718 & & ** &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n0 & FNR & 0.65 & 38 / 8 / 4 & 2.488311e-07 & 2.625475e-07 & & ***** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n0 & FNR & 0.7 & 31 / 11 / 8 & 0.002856943 & 0.00196594 & & ** &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n0 & FNR & 0.75 & 30 / 15 / 5 & 0.009215514 & 0.009579334 & & * &  "
[1] ""
[1] "TQMC-wh-n2 vs. TREE-QMC-n0 & FPR & 0.5 & 33 / 9 / 8 & 0.001517974 & 0.0006522302 & & ** &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n0 & FPR & 0.55 & 35 / 6 / 9 & 1.32592e-05 & 4.598149e-06 & & **** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n0 & FPR & 0.6 & 32 / 13 / 5 & 0.0004205377 & 0.0005079325 & & ** &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n0 & FPR & 0.65 & 40 / 7 / 3 & 1.847749e-07 & 1.467844e-07 & & ***** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n0 & FPR & 0.7 & 32 / 13 / 5 & 0.0008614679 & 0.0009241386 & & ** &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n0 & FPR & 0.75 & 36 / 13 / 1 & 0.0005345105 & 0.0005136653 & & ** &  "
[1] ""
[1] ""
[1] "done"
