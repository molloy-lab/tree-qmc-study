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

        writeme <- paste("TQMC-wh-n2 vs.", mthd, "&",
                         emet, "&",
                         as.character(nbps), "&",
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
[1] "TQMC-wh-n2 vs. ASTEROID & FNR & 50 & 13 / 28 / 9 & 0.02229238 & 0.0179618 & & * &  ( ASTEROID better)"
[1] "TQMC-wh-n2 vs. ASTEROID & FNR & 100 & 12 / 27 / 11 & 0.01183831 & 0.01055658 & & * &  ( ASTEROID better)"
[1] "TQMC-wh-n2 vs. ASTEROID & FNR & 200 & 13 / 24 / 13 & 0.0166201 & 0.02508322 & & * &  ( ASTEROID better)"
[1] "TQMC-wh-n2 vs. ASTEROID & FNR & 500 & 15 / 20 / 15 & 0.2314291 & 0.2800246 & &  &  "
[1] ""
[1] "TQMC-wh-n2 vs. ASTEROID & FPR & 50 & 19 / 27 / 4 & 0.2690406 & 0.257953 & &  &  "
[1] "TQMC-wh-n2 vs. ASTEROID & FPR & 100 & 24 / 23 / 3 & 0.6690783 & 0.7071011 & &  &  "
[1] "TQMC-wh-n2 vs. ASTEROID & FPR & 200 & 27 / 19 / 4 & 0.1983391 & 0.1966187 & &  &  "
[1] "TQMC-wh-n2 vs. ASTEROID & FPR & 500 & 32 / 11 / 7 & 0.0008255697 & 0.0006741725 & & ** &  "
[1] ""
[1] ""
[1] "TQMC-wh-n2 vs. ASTER & FNR & 50 & 29 / 11 / 10 & 0.001355278 & 0.001425032 & & ** &  "
[1] "TQMC-wh-n2 vs. ASTER & FNR & 100 & 24 / 19 / 7 & 0.3086718 & 0.3253 & &  &  "
[1] "TQMC-wh-n2 vs. ASTER & FNR & 200 & 17 / 26 / 7 & 0.3790406 & 0.3133161 & &  &  "
[1] "TQMC-wh-n2 vs. ASTER & FNR & 500 & 21 / 17 / 12 & 0.2529458 & 0.3178422 & &  &  "
[1] ""
[1] "TQMC-wh-n2 vs. ASTER & FPR & 50 & 36 / 10 / 4 & 2.282836e-05 & 2.07947e-05 & & **** & MC "
[1] "TQMC-wh-n2 vs. ASTER & FPR & 100 & 28 / 18 / 4 & 0.006314081 & 0.009292476 & & * &  "
[1] "TQMC-wh-n2 vs. ASTER & FPR & 200 & 29 / 18 / 3 & 0.006777602 & 0.008638316 & & * &  "
[1] "TQMC-wh-n2 vs. ASTER & FPR & 500 & 33 / 14 / 3 & 0.0001358139 & 0.0001752773 & & *** &  "
[1] ""
[1] ""
[1] "TQMC-wh-n2 vs. wASTRID & FNR & 50 & 30 / 17 / 3 & 0.0063477 & 0.007453506 & & * &  "
[1] "TQMC-wh-n2 vs. wASTRID & FNR & 100 & 35 / 10 / 5 & 6.915856e-07 & 1.321502e-06 & & ***** & MC "
[1] "TQMC-wh-n2 vs. wASTRID & FNR & 200 & 41 / 6 / 3 & 9.350335e-08 & 6.366088e-08 & & ***** & MC "
[1] "TQMC-wh-n2 vs. wASTRID & FNR & 500 & 42 / 4 / 4 & 6.234586e-10 & 4.192202e-10 & & ***** & MC "
[1] ""
[1] "TQMC-wh-n2 vs. wASTRID & FPR & 50 & 31 / 17 / 2 & 0.0002473213 & 0.0003477531 & & *** &  "
[1] "TQMC-wh-n2 vs. wASTRID & FPR & 100 & 39 / 9 / 2 & 7.50876e-09 & 1.083473e-08 & & ***** & MC "
[1] "TQMC-wh-n2 vs. wASTRID & FPR & 200 & 45 / 5 / 0 & 1.250307e-10 & 1.250307e-10 & & ***** & MC "
[1] "TQMC-wh-n2 vs. wASTRID & FPR & 500 & 45 / 2 / 3 & 4.334311e-12 & 2.515321e-12 & & ***** & MC "
[1] ""
[1] ""
[1] "TQMC-wh-n2 vs. TREE-QMC-n2-shared & FNR & 50 & 41 / 3 / 6 & 1.045919e-11 & 1.682565e-11 & & ***** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n2-shared & FNR & 100 & 37 / 7 / 6 & 2.857598e-08 & 4.626281e-08 & & ***** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n2-shared & FNR & 200 & 37 / 7 / 6 & 6.555047e-08 & 8.737936e-08 & & ***** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n2-shared & FNR & 500 & 40 / 4 / 6 & 4.407866e-08 & 1.569322e-08 & & ***** & MC "
[1] ""
[1] "TQMC-wh-n2 vs. TREE-QMC-n2-shared & FPR & 50 & 42 / 6 / 2 & 7.368328e-12 & 1.389822e-11 & & ***** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n2-shared & FPR & 100 & 39 / 7 / 4 & 7.394078e-09 & 1.012964e-08 & & ***** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n2-shared & FPR & 200 & 38 / 7 / 5 & 1.47074e-07 & 1.380707e-07 & & ***** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n2-shared & FPR & 500 & 40 / 6 / 4 & 1.428947e-08 & 1.249765e-08 & & ***** & MC "
[1] ""
[1] ""
[1] "TQMC-wh-n2 vs. TREE-QMC-n1 & FNR & 50 & 16 / 17 / 17 & 0.2110743 & 0.5456926 & &  &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n1 & FNR & 100 & 11 / 13 / 26 & 0.903972 & 0.8266549 & &  &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n1 & FNR & 200 & 14 / 14 / 22 & 0.9688128 & 0.9871731 & &  &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n1 & FNR & 500 & 6 / 13 / 31 & 0.0832634 & 0.08807373 & &  &  "
[1] ""
[1] "TQMC-wh-n2 vs. TREE-QMC-n1 & FPR & 50 & 16 / 17 / 17 & 0.4551586 & 0.7416099 & &  &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n1 & FPR & 100 & 13 / 13 / 24 & 0.8950852 & 0.9571749 & &  &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n1 & FPR & 200 & 14 / 14 / 22 & 0.6407185 & 0.8354281 & &  &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n1 & FPR & 500 & 6 / 13 / 31 & 0.01633453 & 0.06314468 & &  &  "
[1] ""
[1] ""
[1] "TQMC-wh-n2 vs. TREE-QMC-n0 & FNR & 50 & 31 / 9 / 10 & 0.0001431637 & 0.0001300513 & & *** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n0 & FNR & 100 & 30 / 12 / 8 & 0.0003918257 & 0.0005882718 & & ** &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n0 & FNR & 200 & 26 / 11 / 13 & 0.003317434 & 0.004399326 & & ** &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n0 & FNR & 500 & 27 / 8 / 15 & 0.001571823 & 0.0009226777 & & ** &  "
[1] ""
[1] "TQMC-wh-n2 vs. TREE-QMC-n0 & FPR & 50 & 32 / 9 / 9 & 4.512263e-05 & 4.787886e-05 & & **** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n0 & FPR & 100 & 32 / 13 / 5 & 0.0004205377 & 0.0005079325 & & ** &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n0 & FPR & 200 & 28 / 13 / 9 & 0.003030112 & 0.004133546 & & ** &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n0 & FPR & 500 & 28 / 9 / 13 & 0.0003299674 & 0.0004092083 & & *** &  "
[1] ""
[1] ""
[1] "done"