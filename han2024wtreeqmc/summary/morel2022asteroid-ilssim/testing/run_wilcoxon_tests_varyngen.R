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

        writeme <- paste("TQMC-wh-n2 vs.", mthd, "&",
                         emet, "&",
                         as.character(ngen), "&",
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
[1] "TQMC-wh-n2 vs. ASTEROID & FNR & 250 & 12 / 36 / 2 & 2.481273e-06 & 3.180854e-06 & & ***** & MC ( ASTEROID better)"
[1] "TQMC-wh-n2 vs. ASTEROID & FNR & 500 & 15 / 28 / 7 & 0.002102961 & 0.003828139 & & ** &  ( ASTEROID better)"
[1] "TQMC-wh-n2 vs. ASTEROID & FNR & 1000 & 12 / 27 / 11 & 0.01183831 & 0.01055658 & & * &  ( ASTEROID better)"
[1] "TQMC-wh-n2 vs. ASTEROID & FNR & 2000 & 13 / 15 / 22 & 0.7740243 & 0.728005 & &  &  "
[1] ""
[1] "TQMC-wh-n2 vs. ASTEROID & FPR & 250 & 28 / 22 / 0 & 0.8538695 & 0.8538695 & &  &  "
[1] "TQMC-wh-n2 vs. ASTEROID & FPR & 500 & 26 / 23 / 1 & 0.2956459 & 0.3033265 & &  &  "
[1] "TQMC-wh-n2 vs. ASTEROID & FPR & 1000 & 24 / 23 / 3 & 0.6690783 & 0.7071011 & &  &  "
[1] "TQMC-wh-n2 vs. ASTEROID & FPR & 2000 & 21 / 13 / 16 & 0.05640108 & 0.08508092 & &  &  "
[1] ""
[1] ""
[1] "TQMC-wh-n2 vs. ASTER & FNR & 250 & 28 / 16 / 6 & 0.1366722 & 0.1138653 & &  &  "
[1] "TQMC-wh-n2 vs. ASTER & FNR & 500 & 30 / 14 / 6 & 0.005711896 & 0.005891858 & & * &  "
[1] "TQMC-wh-n2 vs. ASTER & FNR & 1000 & 24 / 19 / 7 & 0.3086718 & 0.3253 & &  &  "
[1] "TQMC-wh-n2 vs. ASTER & FNR & 2000 & 34 / 9 / 7 & 0.001270064 & 0.0005541056 & & ** &  "
[1] ""
[1] "TQMC-wh-n2 vs. ASTER & FPR & 250 & 40 / 10 / 0 & 2.245184e-08 & 2.245184e-08 & & ***** & MC "
[1] "TQMC-wh-n2 vs. ASTER & FPR & 500 & 39 / 11 / 0 & 2.988292e-07 & 2.988292e-07 & & ***** & MC "
[1] "TQMC-wh-n2 vs. ASTER & FPR & 1000 & 28 / 18 / 4 & 0.006314081 & 0.009292476 & & * &  "
[1] "TQMC-wh-n2 vs. ASTER & FPR & 2000 & 39 / 8 / 3 & 1.796394e-06 & 1.407976e-06 & & ***** & MC "
[1] ""
[1] ""
[1] "TQMC-wh-n2 vs. wASTRID & FNR & 250 & 38 / 9 / 3 & 1.261953e-05 & 9.912823e-06 & & **** & MC "
[1] "TQMC-wh-n2 vs. wASTRID & FNR & 500 & 40 / 7 / 3 & 6.420677e-08 & 5.699137e-08 & & ***** & MC "
[1] "TQMC-wh-n2 vs. wASTRID & FNR & 1000 & 35 / 10 / 5 & 6.915856e-07 & 1.321502e-06 & & ***** & MC "
[1] "TQMC-wh-n2 vs. wASTRID & FNR & 2000 & 36 / 8 / 6 & 2.295702e-07 & 3.496234e-07 & & ***** & MC "
[1] ""
[1] "TQMC-wh-n2 vs. wASTRID & FPR & 250 & 45 / 5 / 0 & 1.44329e-10 & 1.44329e-10 & & ***** & MC "
[1] "TQMC-wh-n2 vs. wASTRID & FPR & 500 & 43 / 7 / 0 & 1.427356e-10 & 1.427356e-10 & & ***** & MC "
[1] "TQMC-wh-n2 vs. wASTRID & FPR & 1000 & 39 / 9 / 2 & 7.50876e-09 & 1.083473e-08 & & ***** & MC "
[1] "TQMC-wh-n2 vs. wASTRID & FPR & 2000 & 40 / 6 / 4 & 9.542873e-10 & 1.343324e-09 & & ***** & MC "
[1] ""
[1] ""
[1] "TQMC-wh-n2 vs. TREE-QMC-n2-shared & FNR & 250 & 41 / 5 / 4 & 3.554419e-10 & 4.049241e-10 & & ***** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n2-shared & FNR & 500 & 41 / 6 / 3 & 2.374492e-10 & 3.431069e-10 & & ***** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n2-shared & FNR & 1000 & 37 / 7 / 6 & 2.857598e-08 & 4.626281e-08 & & ***** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n2-shared & FNR & 2000 & 41 / 4 / 5 & 5.780976e-11 & 7.747758e-11 & & ***** & MC "
[1] ""
[1] "TQMC-wh-n2 vs. TREE-QMC-n2-shared & FPR & 250 & 41 / 8 / 1 & 6.221974e-10 & 7.553567e-10 & & ***** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n2-shared & FPR & 500 & 45 / 5 / 0 & 8.795276e-11 & 8.795276e-11 & & ***** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n2-shared & FPR & 1000 & 39 / 7 / 4 & 7.394078e-09 & 1.012964e-08 & & ***** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n2-shared & FPR & 2000 & 42 / 4 / 4 & 2.54488e-10 & 1.98412e-10 & & ***** & MC "
[1] ""
[1] ""
[1] "TQMC-wh-n2 vs. TREE-QMC-n1 & FNR & 250 & 18 / 14 / 18 & 0.5322092 & 0.4932696 & &  &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n1 & FNR & 500 & 20 / 13 / 17 & 0.1522265 & 0.1693771 & &  &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n1 & FNR & 1000 & 11 / 13 / 26 & 0.903972 & 0.8266549 & &  &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n1 & FNR & 2000 & 14 / 8 / 28 & 0.1527262 & 0.1716685 & &  &  "
[1] ""
[1] "TQMC-wh-n2 vs. TREE-QMC-n1 & FPR & 250 & 16 / 19 / 15 & 0.737086 & 0.9903226 & &  &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n1 & FPR & 500 & 22 / 16 / 12 & 0.1925415 & 0.2214188 & &  &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n1 & FPR & 1000 & 13 / 13 / 24 & 0.8950852 & 0.9571749 & &  &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n1 & FPR & 2000 & 14 / 8 / 28 & 0.1782484 & 0.1807036 & &  &  "
[1] ""
[1] ""
[1] "TQMC-wh-n2 vs. TREE-QMC-n0 & FNR & 250 & 34 / 8 / 8 & 1.485849e-06 & 2.166327e-06 & & ***** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n0 & FNR & 500 & 33 / 10 / 7 & 6.242867e-06 & 1.142868e-05 & & **** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n0 & FNR & 1000 & 30 / 12 / 8 & 0.0003918257 & 0.0005882718 & & ** &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n0 & FNR & 2000 & 36 / 8 / 6 & 1.561854e-05 & 1.001715e-05 & & **** & MC "
[1] ""
[1] "TQMC-wh-n2 vs. TREE-QMC-n0 & FPR & 250 & 35 / 11 / 4 & 1.788362e-06 & 2.951376e-06 & & ***** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n0 & FPR & 500 & 36 / 9 / 5 & 3.574006e-06 & 3.664306e-06 & & ***** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n0 & FPR & 1000 & 32 / 13 / 5 & 0.0004205377 & 0.0005079325 & & ** &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n0 & FPR & 2000 & 36 / 8 / 6 & 1.364693e-05 & 8.909193e-06 & & **** & MC "
[1] ""
[1] ""
[1] "done"
