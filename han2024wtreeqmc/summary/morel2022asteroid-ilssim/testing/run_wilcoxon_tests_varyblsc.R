#!/usr/bin/env Rscript

library(exactRankTests)
library(coin)

data <- read.csv("../csvs/data-varyblsc-for-testing.csv")
data$BLSC <- as.factor(data$BLSC)
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
        for (blsc in blscs) {
        
        xdf <- data[data$BLSC == blsc, ]

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
                         as.character(blsc), "&",
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
[1] "TQMC-wh-n2 vs. ASTEROID & FNR & 0.05 & 14 / 26 / 10 & 0.07256296 & 0.06029934 & &  &  "
[1] "TQMC-wh-n2 vs. ASTEROID & FNR & 0.1 & 9 / 34 / 7 & 0.0006877598 & 0.0003294591 & & ** &  ( ASTEROID better)"
[1] "TQMC-wh-n2 vs. ASTEROID & FNR & 1 & 12 / 27 / 11 & 0.01183831 & 0.01055658 & & * &  ( ASTEROID better)"
[1] "TQMC-wh-n2 vs. ASTEROID & FNR & 10 & 18 / 18 / 14 & 0.5518875 & 0.7110188 & &  &  "
[1] "TQMC-wh-n2 vs. ASTEROID & FNR & 100 & 21 / 25 / 4 & 0.2413859 & 0.2660246 & &  &  "
[1] "TQMC-wh-n2 vs. ASTEROID & FNR & 200 & 12 / 29 / 9 & 0.001171435 & 0.001540183 & & ** &  ( ASTEROID better)"
[1] ""
[1] "TQMC-wh-n2 vs. ASTEROID & FPR & 0.05 & 26 / 23 / 1 & 0.8456339 & 0.8613803 & &  &  "
[1] "TQMC-wh-n2 vs. ASTEROID & FPR & 0.1 & 16 / 29 / 5 & 0.04343865 & 0.04104225 & & * &  ( ASTEROID better)"
[1] "TQMC-wh-n2 vs. ASTEROID & FPR & 1 & 24 / 23 / 3 & 0.6690783 & 0.7071011 & &  &  "
[1] "TQMC-wh-n2 vs. ASTEROID & FPR & 10 & 31 / 16 / 3 & 0.02371947 & 0.02249649 & & * &  "
[1] "TQMC-wh-n2 vs. ASTEROID & FPR & 100 & 25 / 23 / 2 & 0.7086876 & 0.7394276 & &  &  "
[1] "TQMC-wh-n2 vs. ASTEROID & FPR & 200 & 17 / 29 / 4 & 0.01908957 & 0.02131339 & & * &  ( ASTEROID better)"
[1] ""
[1] ""
[1] "TQMC-wh-n2 vs. ASTER & FNR & 0.05 & 27 / 17 / 6 & 0.09089404 & 0.0911233 & &  &  "
[1] "TQMC-wh-n2 vs. ASTER & FNR & 0.1 & 27 / 15 / 8 & 0.05273338 & 0.0495209 & &  &  "
[1] "TQMC-wh-n2 vs. ASTER & FNR & 1 & 24 / 19 / 7 & 0.3086718 & 0.3253 & &  &  "
[1] "TQMC-wh-n2 vs. ASTER & FNR & 10 & 29 / 13 / 8 & 0.0002666164 & 0.0006134766 & & ** &  "
[1] "TQMC-wh-n2 vs. ASTER & FNR & 100 & 37 / 10 / 3 & 1.435091e-06 & 1.671657e-06 & & ***** & MC "
[1] "TQMC-wh-n2 vs. ASTER & FNR & 200 & 25 / 21 / 4 & 0.6237985 & 0.6097468 & &  &  "
[1] ""
[1] "TQMC-wh-n2 vs. ASTER & FPR & 0.05 & 32 / 16 / 2 & 0.002214891 & 0.002416919 & & ** &  "
[1] "TQMC-wh-n2 vs. ASTER & FPR & 0.1 & 29 / 15 / 6 & 0.00114547 & 0.002004851 & & ** &  "
[1] "TQMC-wh-n2 vs. ASTER & FPR & 1 & 28 / 18 / 4 & 0.006314081 & 0.009292476 & & * &  "
[1] "TQMC-wh-n2 vs. ASTER & FPR & 10 & 37 / 12 / 1 & 4.597739e-07 & 5.473287e-07 & & ***** & MC "
[1] "TQMC-wh-n2 vs. ASTER & FPR & 100 & 41 / 9 / 0 & 7.629392e-08 & 7.629392e-08 & & ***** & MC "
[1] "TQMC-wh-n2 vs. ASTER & FPR & 200 & 28 / 20 / 2 & 0.1815662 & 0.1820177 & &  &  "
[1] ""
[1] ""
[1] "TQMC-wh-n2 vs. wASTRID & FNR & 0.05 & 30 / 14 / 6 & 0.005101829 & 0.005369389 & & * &  "
[1] "TQMC-wh-n2 vs. wASTRID & FNR & 0.1 & 35 / 10 / 5 & 8.425792e-06 & 9.914917e-06 & & **** & MC "
[1] "TQMC-wh-n2 vs. wASTRID & FNR & 1 & 35 / 10 / 5 & 6.915856e-07 & 1.321502e-06 & & ***** & MC "
[1] "TQMC-wh-n2 vs. wASTRID & FNR & 10 & 34 / 13 / 3 & 1.793834e-05 & 2.565705e-05 & & **** & MC "
[1] "TQMC-wh-n2 vs. wASTRID & FNR & 100 & 34 / 14 / 2 & 0.0006888812 & 0.0007022983 & & ** &  "
[1] "TQMC-wh-n2 vs. wASTRID & FNR & 200 & 26 / 23 / 1 & 0.585614 & 0.5867599 & &  &  "
[1] ""
[1] "TQMC-wh-n2 vs. wASTRID & FPR & 0.05 & 35 / 14 / 1 & 0.0003081293 & 0.0003138559 & & *** &  "
[1] "TQMC-wh-n2 vs. wASTRID & FPR & 0.1 & 38 / 10 / 2 & 1.230211e-07 & 1.581186e-07 & & ***** & MC "
[1] "TQMC-wh-n2 vs. wASTRID & FPR & 1 & 39 / 9 / 2 & 7.50876e-09 & 1.083473e-08 & & ***** & MC "
[1] "TQMC-wh-n2 vs. wASTRID & FPR & 10 & 39 / 10 / 1 & 9.903916e-08 & 1.101679e-07 & & ***** & MC "
[1] "TQMC-wh-n2 vs. wASTRID & FPR & 100 & 35 / 14 / 1 & 0.0002072638 & 0.0002140064 & & *** &  "
[1] "TQMC-wh-n2 vs. wASTRID & FPR & 200 & 27 / 23 / 0 & 0.2581312 & 0.2581312 & &  &  "
[1] ""
[1] ""
[1] "TQMC-wh-n2 vs. TREE-QMC-n2-shared & FNR & 0.05 & 49 / 1 / 0 & 5.329071e-15 & 5.329071e-15 & & ***** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n2-shared & FNR & 0.1 & 48 / 2 / 0 & 9.947598e-14 & 9.947598e-14 & & ***** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n2-shared & FNR & 1 & 37 / 7 / 6 & 2.857598e-08 & 4.626281e-08 & & ***** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n2-shared & FNR & 10 & 43 / 6 / 1 & 2.420144e-10 & 2.510667e-10 & & ***** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n2-shared & FNR & 100 & 48 / 2 / 0 & 2.131628e-14 & 2.131628e-14 & & ***** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n2-shared & FNR & 200 & 46 / 4 / 0 & 9.84457e-12 & 9.84457e-12 & & ***** & MC "
[1] ""
[1] "TQMC-wh-n2 vs. TREE-QMC-n2-shared & FPR & 0.05 & 49 / 1 / 0 & 3.552714e-15 & 3.552714e-15 & & ***** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n2-shared & FPR & 0.1 & 48 / 2 / 0 & 5.861978e-14 & 5.861978e-14 & & ***** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n2-shared & FPR & 1 & 39 / 7 / 4 & 7.394078e-09 & 1.012964e-08 & & ***** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n2-shared & FPR & 10 & 44 / 6 / 0 & 5.68825e-11 & 5.68825e-11 & & ***** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n2-shared & FPR & 100 & 49 / 1 / 0 & 1.24345e-14 & 1.24345e-14 & & ***** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n2-shared & FPR & 200 & 47 / 3 / 0 & 9.082513e-12 & 9.082513e-12 & & ***** & MC "
[1] ""
[1] ""
[1] "TQMC-wh-n2 vs. TREE-QMC-n1 & FNR & 0.05 & 15 / 24 / 11 & 0.3062118 & 0.2335955 & &  &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n1 & FNR & 0.1 & 14 / 22 / 14 & 0.1544441 & 0.1523943 & &  &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n1 & FNR & 1 & 11 / 13 / 26 & 0.903972 & 0.8266549 & &  &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n1 & FNR & 10 & 21 / 10 / 19 & 0.005176618 & 0.0144605 & & * &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n1 & FNR & 100 & 30 / 12 / 8 & 0.002448538 & 0.002319233 & & ** &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n1 & FNR & 200 & 22 / 19 / 9 & 0.6496359 & 0.6387548 & &  &  "
[1] ""
[1] "TQMC-wh-n2 vs. TREE-QMC-n1 & FPR & 0.05 & 16 / 25 / 9 & 0.352719 & 0.2773672 & &  &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n1 & FPR & 0.1 & 14 / 23 / 13 & 0.1849054 & 0.1537121 & &  &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n1 & FPR & 1 & 13 / 13 / 24 & 0.8950852 & 0.9571749 & &  &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n1 & FPR & 10 & 24 / 10 / 16 & 0.02129918 & 0.01496583 & & * &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n1 & FPR & 100 & 31 / 12 / 7 & 0.01124468 & 0.007305996 & & * &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n1 & FPR & 200 & 22 / 19 / 9 & 0.4388212 & 0.4782963 & &  &  "
[1] ""
[1] ""
[1] "TQMC-wh-n2 vs. TREE-QMC-n0 & FNR & 0.05 & 29 / 14 / 7 & 0.009422605 & 0.009666751 & & * &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n0 & FNR & 0.1 & 34 / 10 / 6 & 7.137124e-05 & 6.409277e-05 & & *** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n0 & FNR & 1 & 30 / 12 / 8 & 0.0003918257 & 0.0005882718 & & ** &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n0 & FNR & 10 & 34 / 12 / 4 & 8.107025e-06 & 1.311101e-05 & & **** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n0 & FNR & 100 & 37 / 7 / 6 & 1.59838e-08 & 3.025548e-08 & & ***** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n0 & FNR & 200 & 31 / 13 / 6 & 0.0004322548 & 0.0005982124 & & ** &  "
[1] ""
[1] "TQMC-wh-n2 vs. TREE-QMC-n0 & FPR & 0.05 & 33 / 15 / 2 & 0.004233391 & 0.004140666 & & ** &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n0 & FPR & 0.1 & 35 / 12 / 3 & 2.21552e-05 & 2.642192e-05 & & **** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n0 & FPR & 1 & 32 / 13 / 5 & 0.0004205377 & 0.0005079325 & & ** &  "
[1] "TQMC-wh-n2 vs. TREE-QMC-n0 & FPR & 10 & 37 / 12 / 1 & 4.369193e-06 & 4.714765e-06 & & ***** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n0 & FPR & 100 & 40 / 8 / 2 & 2.666075e-08 & 2.896088e-08 & & ***** & MC "
[1] "TQMC-wh-n2 vs. TREE-QMC-n0 & FPR & 200 & 35 / 13 / 2 & 0.0003031529 & 0.0003003401 & & *** &  "
[1] ""
[1] ""
[1] "done"
