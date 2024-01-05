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

        writeme <- paste("TQMC-wh-n2 vs.", mthd, "&",
                         emet, "&",
                         as.character(ntax), "&",
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
[1] "TQMC-wh-n2 vs. ASTEROID & FNR & 25 & 13 / 18 / 19 & 0.2284313 & 0.279966 & &  &  "
[1] "TQMC-wh-n2 vs. ASTEROID & FNR & 75 & 14 / 27 / 9 & 0.004961303 & 0.007824902 & & * &  ( ASTEROID better)"
[1] "TQMC-wh-n2 vs. ASTEROID & FNR & 50 & 12 / 27 / 11 & 0.01183831 & 0.01055658 & & * &  ( ASTEROID better)"
[1] "TQMC-wh-n2 vs. ASTEROID & FNR & 100 & 13 / 28 / 9 & 0.1062621 & 0.06212415 & &  &  "
[1] "TQMC-wh-n2 vs. ASTEROID & FNR & 125 & 13 / 29 / 8 & 0.002018967 & 0.002639422 & & ** &  ( ASTEROID better)"
[1] "TQMC-wh-n2 vs. ASTEROID & FNR & 150 & 15 / 32 / 3 & 0.001363031 & 0.001541501 & & ** &  ( ASTEROID better)"
[1] ""
[1] "TQMC-wh-n2 vs. ASTEROID & FPR & 25 & 20 / 16 / 14 & 0.6948928 & 0.6059339 & &  &  "
[1] "TQMC-wh-n2 vs. ASTEROID & FPR & 75 & 23 / 25 / 2 & 0.739193 & 0.7685263 & &  &  "
[1] "TQMC-wh-n2 vs. ASTEROID & FPR & 50 & 24 / 23 / 3 & 0.6690783 & 0.7071011 & &  &  "
[1] "TQMC-wh-n2 vs. ASTEROID & FPR & 100 & 29 / 20 / 1 & 0.02683605 & 0.02843397 & & * &  "
[1] "TQMC-wh-n2 vs. ASTEROID & FPR & 125 & 30 / 19 / 1 & 0.00905904 & 0.009744651 & & * &  "
[1] "TQMC-wh-n2 vs. ASTEROID & FPR & 150 & 31 / 19 / 0 & 0.107376 & 0.107376 & &  &  "
[1] ""
[1] ""
[1] "TQMC-wh-n2 vs. ASTER & FNR & 25 & 29 / 9 / 12 & 0.000788109 & 0.0005895351 & & ** &  "
[1] "TQMC-wh-n2 vs. ASTER & FNR & 75 & 26 / 18 / 6 & 0.09673945 & 0.1082865 & &  &  "
[1] "TQMC-wh-n2 vs. ASTER & FNR & 50 & 24 / 19 / 7 & 0.3086718 & 0.3253 & &  &  "
[1] "TQMC-wh-n2 vs. ASTER & FNR & 100 & 30 / 15 / 5 & 0.009233083 & 0.009592993 & & * &  "
[1] "TQMC-wh-n2 vs. ASTER & FNR & 125 & 31 / 16 / 3 & 0.002458039 & 0.00290118 & & ** &  "
[1] "TQMC-wh-n2 vs. ASTER & FNR & 150 & 38 / 6 / 6 & 3.05818e-06 & 1.263203e-06 & & ***** & MC "
[1] ""
[1] "TQMC-wh-n2 vs. ASTER & FPR & 25 & 36 / 7 / 7 & 8.042848e-06 & 4.286517e-06 & & **** & MC "
[1] "TQMC-wh-n2 vs. ASTER & FPR & 75 & 33 / 16 / 1 & 0.0002887686 & 0.0003205901 & & *** &  "
[1] "TQMC-wh-n2 vs. ASTER & FPR & 50 & 28 / 18 / 4 & 0.006314081 & 0.009292476 & & * &  "
[1] "TQMC-wh-n2 vs. ASTER & FPR & 100 & 39 / 10 / 1 & 3.727569e-08 & 4.340035e-08 & & ***** & MC "
[1] "TQMC-wh-n2 vs. ASTER & FPR & 125 & 41 / 9 / 0 & 1.887022e-08 & 1.887022e-08 & & ***** & MC "
[1] "TQMC-wh-n2 vs. ASTER & FPR & 150 & 46 / 4 / 0 & 3.677059e-13 & 3.677059e-13 & & ***** & MC "
[1] ""
[1] ""
[1] "TQMC-wh-n2 vs. wASTRID & FNR & 25 & 25 / 14 / 11 & 0.0366626 & 0.04129837 & & * &  "
[1] "TQMC-wh-n2 vs. wASTRID & FNR & 75 & 36 / 9 / 5 & 2.9034e-07 & 4.79056e-07 & & ***** & MC "
[1] "TQMC-wh-n2 vs. wASTRID & FNR & 50 & 35 / 10 / 5 & 6.915856e-07 & 1.321502e-06 & & ***** & MC "
[1] "TQMC-wh-n2 vs. wASTRID & FNR & 100 & 40 / 8 / 2 & 4.746632e-09 & 6.052488e-09 & & ***** & MC "
[1] "TQMC-wh-n2 vs. wASTRID & FNR & 125 & 45 / 4 / 1 & 2.761524e-11 & 2.510703e-11 & & ***** & MC "
[1] "TQMC-wh-n2 vs. wASTRID & FNR & 150 & 48 / 2 / 0 & 4.565237e-13 & 4.565237e-13 & & ***** & MC "
[1] ""
[1] "TQMC-wh-n2 vs. wASTRID & FPR & 25 & 31 / 12 / 7 & 0.002056411 & 0.001819913 & & ** &  "
[1] "TQMC-wh-n2 vs. wASTRID & FPR & 75 & 42 / 7 / 1 & 8.400995e-10 & 9.076189e-10 & & ***** & MC "
[1] "TQMC-wh-n2 vs. wASTRID & FPR & 50 & 39 / 9 / 2 & 7.50876e-09 & 1.083473e-08 & & ***** & MC "
[1] "TQMC-wh-n2 vs. wASTRID & FPR & 100 & 43 / 6 / 1 & 8.277823e-12 & 1.071854e-11 & & ***** & MC "
[1] "TQMC-wh-n2 vs. wASTRID & FPR & 125 & 47 / 3 / 0 & 1.953993e-13 & 1.953993e-13 & & ***** & MC "
[1] "TQMC-wh-n2 vs. wASTRID & FPR & 150 & 49 / 1 / 0 & 8.881784e-15 & 8.881784e-15 & & ***** & MC "
[1] ""
[1] ""
[1] "TQMC-wh-n2 vs. TREE-QMC-n2-shared & FNR & 25 & 28 / 9 / 13 & 2.926651e-05 & 0.0001020106 & & *** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n2-shared & FNR & 75 & 46 / 2 / 2 & 1.350031e-12 & 9.166001e-13 & & ***** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n2-shared & FNR & 50 & 37 / 7 / 6 & 2.857598e-08 & 4.626281e-08 & & ***** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n2-shared & FNR & 100 & 50 / 0 / 0 & 1.776357e-15 & 1.776357e-15 & & ***** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n2-shared & FNR & 125 & 50 / 0 / 0 & 1.776357e-15 & 1.776357e-15 & & ***** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n2-shared & FNR & 150 & 50 / 0 / 0 & 1.776357e-15 & 1.776357e-15 & & ***** & MC "
[1] ""
[1] "TQMC-wh-n2 vs. TREE-QMC-n2-shared & FPR & 25 & 28 / 11 / 11 & 3.397626e-05 & 0.0001708718 & & *** &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n2-shared & FPR & 75 & 46 / 3 / 1 & 7.354117e-13 & 7.212009e-13 & & ***** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n2-shared & FPR & 50 & 39 / 7 / 4 & 7.394078e-09 & 1.012964e-08 & & ***** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n2-shared & FPR & 100 & 50 / 0 / 0 & 1.776357e-15 & 1.776357e-15 & & ***** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n2-shared & FPR & 125 & 50 / 0 / 0 & 1.776357e-15 & 1.776357e-15 & & ***** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n2-shared & FPR & 150 & 50 / 0 / 0 & 1.776357e-15 & 1.776357e-15 & & ***** & MC "
[1] ""
[1] ""
[1] "TQMC-wh-n2 vs. TREE-QMC-n1 & FNR & 25 & 7 / 8 / 35 & 0.9657593 & 0.8944092 & &  &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n1 & FNR & 75 & 22 / 15 / 13 & 0.1283339 & 0.1559355 & &  &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n1 & FNR & 50 & 11 / 13 / 26 & 0.903972 & 0.8266549 & &  &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n1 & FNR & 100 & 21 / 17 / 12 & 0.1570862 & 0.2381887 & &  &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n1 & FNR & 125 & 26 / 16 / 8 & 0.03523107 & 0.04402304 & & * &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n1 & FNR & 150 & 27 / 16 / 7 & 0.02284229 & 0.02837261 & & * &  "
[1] ""
[1] "TQMC-wh-n2 vs. TREE-QMC-n1 & FPR & 25 & 7 / 8 / 35 & 0.9901733 & 0.8679199 & &  &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n1 & FPR & 75 & 22 / 16 / 12 & 0.3217833 & 0.3108119 & &  &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n1 & FPR & 50 & 13 / 13 / 24 & 0.8950852 & 0.9571749 & &  &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n1 & FPR & 100 & 21 / 17 / 12 & 0.1154733 & 0.1984497 & &  &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n1 & FPR & 125 & 27 / 16 / 7 & 0.1582689 & 0.1329282 & &  &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n1 & FPR & 150 & 28 / 16 / 6 & 0.02306327 & 0.02582644 & & * &  "
[1] ""
[1] ""
[1] "TQMC-wh-n2 vs. TREE-QMC-n0 & FNR & 25 & 29 / 8 / 13 & 2.910336e-05 & 5.518943e-05 & & *** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n0 & FNR & 75 & 30 / 11 / 9 & 0.000167518 & 0.000275696 & & *** &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n0 & FNR & 50 & 30 / 12 / 8 & 0.0003918257 & 0.0005882718 & & ** &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n0 & FNR & 100 & 39 / 7 / 4 & 5.016631e-09 & 7.349655e-09 & & ***** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n0 & FNR & 125 & 39 / 6 / 5 & 1.326222e-07 & 8.816784e-08 & & ***** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n0 & FNR & 150 & 46 / 3 / 1 & 6.679102e-13 & 7.460699e-13 & & ***** & MC "
[1] ""
[1] "TQMC-wh-n2 vs. TREE-QMC-n0 & FPR & 25 & 30 / 10 / 10 & 0.0001998023 & 0.0002538633 & & *** &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n0 & FPR & 75 & 30 / 13 / 7 & 0.000218137 & 0.000424311 & & *** &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n0 & FPR & 50 & 32 / 13 / 5 & 0.0004205377 & 0.0005079325 & & ** &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n0 & FPR & 100 & 40 / 8 / 2 & 2.386095e-09 & 3.25587e-09 & & ***** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n0 & FPR & 125 & 42 / 5 / 3 & 5.350343e-08 & 3.024985e-08 & & ***** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n0 & FPR & 150 & 47 / 3 / 0 & 7.638334e-14 & 7.638334e-14 & & ***** & MC "
[1] ""
[1] ""
[1] "done"

