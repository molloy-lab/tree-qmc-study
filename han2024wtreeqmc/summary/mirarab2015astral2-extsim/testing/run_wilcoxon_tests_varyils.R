#!/usr/bin/env Rscript

library(exactRankTests)
library(coin)

data <- read.csv("../csvs/data-varyils-for-testing.csv")

data$NTAX <- as.factor(data$NTAX)
data$HGHT <- as.factor(data$HGHT)
data$RATE <- as.factor(data$RATE)
data$SUPP <- as.factor(data$SUPP)
data$NGEN <- as.factor(data$NGEN)
data$REPL <- as.factor(data$REPL)

#data <- subset(data, NTAX == 200)

supps <- c("abayes")
ngens <- c(50, 200, 1000)

hghts <- c("500000", "2000000", "10000000")
ilss <- c("0.25X", "1X", "5X")

rates <- c("1e-07", "1e-06")
specs <- c("deep", "shallow")

mthds <- c("CAML", "ASTER-h", "wASTRID")

ntests <- length(supps) * length(ngens) * 11 * length(mthds)
threshold0 <- 0.05
threshold_bonferoni <- threshold0 / ntests

writeme <- paste("alpha =", as.character(threshold0))
print(writeme)
print(paste("bonferroni alpha =", as.character(threshold_bonferoni), 
            "=", as.character(threshold0), "/", as.character(ntests)))

for (mthd in mthds) {
    for (supp in supps) {
        for (i in c(1, 2, 3)) {
                hght <- hghts[i]
                ils <- ilss[i]
            for (j in c(1, 2)) {
                rate <- rates[j]
                spec <- specs[j]

                for (ngen in ngens) {

        #print(paste(as.character(i), as.character(j)))

        df <- data[data$SUPP == supp, ]
        df <- df[df$NGEN == ngen, ]
        df <- df[df$HGHT == hght, ]
        df <- df[df$RATE == rate, ]

        # Find ties
        mthd1 <- df$TQMCwhn2xSEFN
        if (mthd == "CAML") {
            mthd2 <- df$CAMLxSEFN
        } else if (mthd == "ASTER-h") {
            mthd2 <- df$ASTERHxSEFN
        } else if (mthd == "wASTRID") {
            mthd2 <- df$WASTRIDxSEFN
        } else {
            print("ERROR")
            exit()
        }
        diff <- mthd2 - mthd1  # positive means TQMC-wh-n2 is better

        # Count ties
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
                         supp, "&",
                         ils, "&",
                         spec, "&",
                         as.character(ngen), "&",
                         as.character(ntreeqmc), "/",
                         as.character(nother), "/",
                         as.character(ntie), "&",
                         format(xwsr, scientific = TRUE), "&", 
                         format(ywsr, scientific = TRUE), "&",
                         "&", stars, "&", mc, note)
        print(writeme)

                }
            }
        }
    }
}

print("done")
exit()

[1] "alpha = 0.05"
[1] "bonferroni alpha = 0.000505050505050505 = 0.05 / 99"
[1] "TQMC-wh-n2 vs. CAML & abayes & 0.25X & deep & 50 & 50 / 0 / 0 & 1.776357e-15 & 1.776357e-15 & & ***** & MC "
[1] "TQMC-wh-n2 vs. CAML & abayes & 0.25X & deep & 200 & 50 / 0 / 0 & 1.776357e-15 & 1.776357e-15 & & ***** & MC "
[1] "TQMC-wh-n2 vs. CAML & abayes & 0.25X & deep & 1000 & 46 / 4 / 0 & 7.654433e-08 & 7.654433e-08 & & ***** & MC "
[1] "TQMC-wh-n2 vs. CAML & abayes & 0.25X & shallow & 50 & 47 / 0 / 0 & 1.421085e-14 & 1.421085e-14 & & ***** & MC "
[1] "TQMC-wh-n2 vs. CAML & abayes & 0.25X & shallow & 200 & 47 / 0 / 0 & 1.421085e-14 & 1.421085e-14 & & ***** & MC "
[1] "TQMC-wh-n2 vs. CAML & abayes & 0.25X & shallow & 1000 & 46 / 0 / 1 & 2.842171e-14 & 2.842171e-14 & & ***** & MC "
[1] "TQMC-wh-n2 vs. CAML & abayes & 1X & deep & 50 & 38 / 10 / 2 & 4.729318e-06 & 4.575677e-06 & & ***** & MC "
[1] "TQMC-wh-n2 vs. CAML & abayes & 1X & deep & 200 & 37 / 9 / 4 & 0.0004051359 & 0.0002279103 & & *** & MC "
[1] "TQMC-wh-n2 vs. CAML & abayes & 1X & deep & 1000 & 24 / 19 / 7 & 0.9784243 & 0.8870725 & &  &  "
[1] "TQMC-wh-n2 vs. CAML & abayes & 1X & shallow & 50 & 44 / 3 / 3 & 3.012701e-12 & 3.538503e-12 & & ***** & MC "
[1] "TQMC-wh-n2 vs. CAML & abayes & 1X & shallow & 200 & 40 / 4 / 6 & 2.048751e-09 & 1.34969e-09 & & ***** & MC "
[1] "TQMC-wh-n2 vs. CAML & abayes & 1X & shallow & 1000 & 33 / 8 / 9 & 0.0003452434 & 0.0001578385 & & *** & MC "
[1] "TQMC-wh-n2 vs. CAML & abayes & 5X & deep & 50 & 12 / 31 / 7 & 0.0009704487 & 0.001014695 & & ** &  ( CAML better)"
[1] "TQMC-wh-n2 vs. CAML & abayes & 5X & deep & 200 & 13 / 28 / 9 & 0.0007806279 & 0.001628425 & & ** &  ( CAML better)"
[1] "TQMC-wh-n2 vs. CAML & abayes & 5X & deep & 1000 & 13 / 29 / 8 & 0.01314316 & 0.01115074 & & * &  ( CAML better)"
[1] "TQMC-wh-n2 vs. CAML & abayes & 5X & shallow & 50 & 36 / 9 / 5 & 5.632589e-06 & 5.51794e-06 & & **** & MC "
[1] "TQMC-wh-n2 vs. CAML & abayes & 5X & shallow & 200 & 38 / 6 / 6 & 8.63588e-06 & 2.94967e-06 & & **** & MC "
[1] "TQMC-wh-n2 vs. CAML & abayes & 5X & shallow & 1000 & 29 / 10 / 11 & 0.005499648 & 0.003169766 & & * &  "
[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 0.25X & deep & 50 & 25 / 21 / 4 & 0.1839068 & 0.2112505 & &  &  "
[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 0.25X & deep & 200 & 24 / 21 / 5 & 0.8416785 & 0.9229223 & &  &  "
[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 0.25X & deep & 1000 & 21 / 17 / 12 & 0.9515619 & 0.7842332 & &  &  "
[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 0.25X & shallow & 50 & 29 / 13 / 5 & 0.01055602 & 0.009546006 & & * &  "
[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 0.25X & shallow & 200 & 29 / 12 / 6 & 0.01537982 & 0.01152872 & & * &  "
[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 0.25X & shallow & 1000 & 26 / 14 / 7 & 0.07923641 & 0.06927969 & &  &  "
[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 1X & deep & 50 & 24 / 14 / 12 & 0.04777693 & 0.05695275 & &  &  "
[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 1X & deep & 200 & 29 / 8 / 13 & 0.0001611809 & 0.0001548738 & & *** & MC "
[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 1X & deep & 1000 & 22 / 11 / 17 & 0.008441326 & 0.0186442 & & * &  "
[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 1X & shallow & 50 & 22 / 17 / 11 & 0.02650681 & 0.07322669 & &  &  "
[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 1X & shallow & 200 & 25 / 10 / 15 & 0.01056484 & 0.00998523 & & * &  "
[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 1X & shallow & 1000 & 17 / 5 / 28 & 0.01054907 & 0.01173639 & & * &  "
[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 5X & deep & 50 & 30 / 3 / 17 & 1.357403e-07 & 2.060551e-07 & & ***** & MC "
[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 5X & deep & 200 & 21 / 11 / 18 & 0.008554659 & 0.02371259 & & * &  "
[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 5X & deep & 1000 & 26 / 2 / 22 & 1.696497e-05 & 2.294779e-06 & & **** & MC "
[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 5X & shallow & 50 & 32 / 10 / 8 & 4.137752e-05 & 6.118376e-05 & & *** & MC "
[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 5X & shallow & 200 & 21 / 13 / 16 & 0.1553411 & 0.1606247 & &  &  "
[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 5X & shallow & 1000 & 17 / 5 / 28 & 0.03724957 & 0.01563883 & & * &  "
[1] "TQMC-wh-n2 vs. wASTRID & abayes & 0.25X & deep & 50 & 35 / 9 / 6 & 9.161746e-06 & 9.006704e-06 & & **** & MC "
[1] "TQMC-wh-n2 vs. wASTRID & abayes & 0.25X & deep & 200 & 37 / 10 / 3 & 1.476171e-05 & 1.351772e-05 & & **** & MC "
[1] "TQMC-wh-n2 vs. wASTRID & abayes & 0.25X & deep & 1000 & 30 / 13 / 7 & 0.1387233 & 0.07745219 & &  &  "
[1] "TQMC-wh-n2 vs. wASTRID & abayes & 0.25X & shallow & 50 & 39 / 5 / 3 & 7.166818e-10 & 8.781171e-10 & & ***** & MC "
[1] "TQMC-wh-n2 vs. wASTRID & abayes & 0.25X & shallow & 200 & 40 / 5 / 2 & 5.238405e-08 & 3.861595e-08 & & ***** & MC "
[1] "TQMC-wh-n2 vs. wASTRID & abayes & 0.25X & shallow & 1000 & 33 / 8 / 6 & 7.542394e-07 & 1.492041e-06 & & ***** & MC "
[1] "TQMC-wh-n2 vs. wASTRID & abayes & 1X & deep & 50 & 25 / 20 / 5 & 0.03693053 & 0.05764566 & &  &  "
[1] "TQMC-wh-n2 vs. wASTRID & abayes & 1X & deep & 200 & 33 / 10 / 7 & 6.022883e-05 & 6.837434e-05 & & *** & MC "
[1] "TQMC-wh-n2 vs. wASTRID & abayes & 1X & deep & 1000 & 30 / 6 / 14 & 4.754454e-05 & 2.847172e-05 & & **** & MC "
[1] "TQMC-wh-n2 vs. wASTRID & abayes & 1X & shallow & 50 & 30 / 13 / 7 & 0.0007217834 & 0.001086662 & & ** &  "
[1] "TQMC-wh-n2 vs. wASTRID & abayes & 1X & shallow & 200 & 26 / 13 / 11 & 0.01845867 & 0.02047486 & & * &  "
[1] "TQMC-wh-n2 vs. wASTRID & abayes & 1X & shallow & 1000 & 26 / 12 / 12 & 0.01553998 & 0.01566447 & & * &  "
[1] "TQMC-wh-n2 vs. wASTRID & abayes & 5X & deep & 50 & 31 / 7 / 12 & 7.58782e-06 & 1.008374e-05 & & **** & MC "
[1] "TQMC-wh-n2 vs. wASTRID & abayes & 5X & deep & 200 & 38 / 5 / 7 & 1.490207e-09 & 2.443812e-09 & & ***** & MC "
[1] "TQMC-wh-n2 vs. wASTRID & abayes & 5X & deep & 1000 & 43 / 2 / 5 & 3.183231e-12 & 3.581135e-12 & & ***** & MC "
[1] "TQMC-wh-n2 vs. wASTRID & abayes & 5X & shallow & 50 & 22 / 13 / 15 & 0.2696331 & 0.1872801 & &  &  "
[1] "TQMC-wh-n2 vs. wASTRID & abayes & 5X & shallow & 200 & 12 / 20 / 18 & 0.2902047 & 0.20998 & &  &  "
[1] "TQMC-wh-n2 vs. wASTRID & abayes & 5X & shallow & 1000 & 14 / 7 / 29 & 0.2642612 & 0.160861 & &  &  "
[1] "done"
