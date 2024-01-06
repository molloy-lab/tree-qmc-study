#!/usr/bin/env Rscript

# install.packages("coin", dependencies = TRUE)

library(exactRankTests)
library(coin)

#https://sphweb.bumc.bu.edu/otlt/mph-modules/bs/bs704_nonparametric/BS704_Nonparametric6.html
#https://stat.ethz.ch/pipermail/r-help/2011-April/274931.html
#https://www.graphpad.com/guides/prism/latest/statistics/stat_pratt.htm
#https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/wilcox.test
# https://www.rdocumentation.org/packages/rstatix/versions/0.7.2/topics/wilcox_test

data <- read.csv("../csvs/data-for-testing.csv")
data$NBPS <- as.factor(data$NBPS)
data$NGEN <- as.factor(data$NGEN)
data$SUPP <- as.factor(data$SUPP)
data$REPL <- as.factor(data$REPL)

ngens <- c("50", "200", "500", "1000")
nbpss <- c("200", "400", "800", "1600")
mthds <- c("ASTER-wh", "ASTRID-ws", "TQMC-n2-refinedpoly")

threshold0 <- 0.05
ntests <- length(ngens) * length(nbpss) * 5
threshold_bonferoni <- threshold0 / ntests

writeme <- paste("alpha =", as.character(threshold0))
print(writeme)
print(paste("bonferroni alpha =", as.character(threshold_bonferoni), 
            "=", as.character(threshold0), "/", as.character(ntests)))

for (mthd in mthds) {

    if (mthd == "TQMC-n2-refinedpoly") {
        supps <- c("abayes")  # abayes trees have polys refined
    } else {
        supps <- c("abayes", "bs")
    }

    for (supp in supps) {
        for (ngen in ngens) {
            for (nbps in nbpss) {
        
        df <- data[data$SUPP == supp, ]
        df <- df[df$NGEN == ngen, ]
        df <- df[df$NBPS == nbps, ]

        # Set-up
        mthd1 <- df$TQMCwhn2xSERF #FN
        if (mthd == "ASTER-wh") {
            mthd2 <- df$ASTERHxSERF #FN
            supplab = supp
        } else if (mthd == "ASTRID-ws") {
            mthd2 <- df$WASTRIDxSERF #FN
            supplab = supp
        } else if (mthd == "TQMC-n2-refinedpoly") {
            mthd2 <- df$TQMCwnn2xSERF #FN
            supplab = "none"
        } else {
            print("ERROR")
            exit()
        }

        # Count ties
        diff <- mthd2 - mthd1  # lower is better
        ntie <- sum(diff == 0)
        nother <- sum(diff < 0)
        ntreeqmc <- sum(diff > 0)

        # Perform wilcoxon test
        #xwsr <-  wilcox.exact(mthd1, mthd2,
        #                      paired = TRUE,
        #                      alternative = "two.sided",
        #                      conf.int = TRUE)
        #xwsr <- xwsr$p.value
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

        writeme <- paste("TQMC-wh_n2 vs.", mthd, "&",
                         supplab, "& &",
                         as.character(ngen), "&",
                         as.character(nbps), "& &",
                         as.character(ntreeqmc), "&",
                         as.character(nother), "&",
                         as.character(ntie), "& &",
                         format(wsr, scientific = TRUE), "&", 
                         stars, "&", mc, "&", note, " \\\\")
        message(writeme)

            }
        }
    }
    message("")
}

print("done")
exit()

[1] "alpha = 0.05"
[1] "bonferroni alpha = 0.000625 = 0.05 / 80"
TQMC-wh_n2 vs. ASTER-wh & abayes & & 50 & 200 & & 23 & 17 & 10 & & 0.3341381 &  &  &   \\
TQMC-wh_n2 vs. ASTER-wh & abayes & & 50 & 400 & & 33 & 7 & 10 & & 0.000387488 & *** & MC &   \\
TQMC-wh_n2 vs. ASTER-wh & abayes & & 50 & 800 & & 27 & 8 & 15 & & 0.0002738126 & *** & MC &   \\
TQMC-wh_n2 vs. ASTER-wh & abayes & & 50 & 1600 & & 23 & 6 & 21 & & 0.0004458502 & *** & MC &   \\
TQMC-wh_n2 vs. ASTER-wh & abayes & & 200 & 200 & & 30 & 9 & 11 & & 0.006861866 & * &  &   \\
TQMC-wh_n2 vs. ASTER-wh & abayes & & 200 & 400 & & 21 & 14 & 15 & & 0.1534584 &  &  &   \\
TQMC-wh_n2 vs. ASTER-wh & abayes & & 200 & 800 & & 18 & 11 & 21 & & 0.1550322 &  &  &   \\
TQMC-wh_n2 vs. ASTER-wh & abayes & & 200 & 1600 & & 22 & 10 & 18 & & 0.02215265 & * &  &   \\
TQMC-wh_n2 vs. ASTER-wh & abayes & & 500 & 200 & & 22 & 15 & 13 & & 0.7603878 &  &  &   \\
TQMC-wh_n2 vs. ASTER-wh & abayes & & 500 & 400 & & 14 & 22 & 14 & & 0.529488 &  &  &   \\
TQMC-wh_n2 vs. ASTER-wh & abayes & & 500 & 800 & & 18 & 11 & 21 & & 0.5209156 &  &  &   \\
TQMC-wh_n2 vs. ASTER-wh & abayes & & 500 & 1600 & & 21 & 13 & 16 & & 0.1982383 &  &  &   \\
TQMC-wh_n2 vs. ASTER-wh & abayes & & 1000 & 200 & & 20 & 9 & 21 & & 0.006101582 & * &  &   \\
TQMC-wh_n2 vs. ASTER-wh & abayes & & 1000 & 400 & & 20 & 6 & 24 & & 0.003370821 & ** &  &   \\
TQMC-wh_n2 vs. ASTER-wh & abayes & & 1000 & 800 & & 20 & 7 & 23 & & 0.007649198 & * &  &   \\
TQMC-wh_n2 vs. ASTER-wh & abayes & & 1000 & 1600 & & 19 & 5 & 26 & & 0.001789212 & ** &  &   \\
TQMC-wh_n2 vs. ASTER-wh & bs & & 50 & 200 & & 22 & 17 & 11 & & 0.3525142 &  &  &   \\
TQMC-wh_n2 vs. ASTER-wh & bs & & 50 & 400 & & 19 & 14 & 17 & & 0.4983073 &  &  &   \\
TQMC-wh_n2 vs. ASTER-wh & bs & & 50 & 800 & & 23 & 18 & 9 & & 0.8386277 &  &  &   \\
TQMC-wh_n2 vs. ASTER-wh & bs & & 50 & 1600 & & 21 & 16 & 13 & & 0.1983822 &  &  &   \\
TQMC-wh_n2 vs. ASTER-wh & bs & & 200 & 200 & & 19 & 17 & 14 & & 0.6023837 &  &  &   \\
TQMC-wh_n2 vs. ASTER-wh & bs & & 200 & 400 & & 21 & 11 & 18 & & 0.06322916 &  &  &   \\
TQMC-wh_n2 vs. ASTER-wh & bs & & 200 & 800 & & 17 & 18 & 15 & & 0.6384818 &  &  &   \\
TQMC-wh_n2 vs. ASTER-wh & bs & & 200 & 1600 & & 17 & 11 & 22 & & 0.2262515 &  &  &   \\
TQMC-wh_n2 vs. ASTER-wh & bs & & 500 & 200 & & 23 & 18 & 9 & & 0.7983978 &  &  &   \\
TQMC-wh_n2 vs. ASTER-wh & bs & & 500 & 400 & & 21 & 14 & 15 & & 0.112634 &  &  &   \\
TQMC-wh_n2 vs. ASTER-wh & bs & & 500 & 800 & & 17 & 10 & 23 & & 0.3459887 &  &  &   \\
TQMC-wh_n2 vs. ASTER-wh & bs & & 500 & 1600 & & 13 & 8 & 29 & & 0.2260199 &  &  &   \\
TQMC-wh_n2 vs. ASTER-wh & bs & & 1000 & 200 & & 16 & 19 & 15 & & 0.6453902 &  &  &   \\
TQMC-wh_n2 vs. ASTER-wh & bs & & 1000 & 400 & & 16 & 10 & 24 & & 0.1127932 &  &  &   \\
TQMC-wh_n2 vs. ASTER-wh & bs & & 1000 & 800 & & 17 & 6 & 27 & & 0.08854365 &  &  &   \\
TQMC-wh_n2 vs. ASTER-wh & bs & & 1000 & 1600 & & 13 & 10 & 27 & & 0.4147711 &  &  &   \\

TQMC-wh_n2 vs. ASTRID-ws & abayes & & 50 & 200 & & 30 & 16 & 4 & & 0.02684553 & * &  &   \\
TQMC-wh_n2 vs. ASTRID-ws & abayes & & 50 & 400 & & 33 & 12 & 5 & & 0.001025001 & ** &  &   \\
TQMC-wh_n2 vs. ASTRID-ws & abayes & & 50 & 800 & & 35 & 10 & 5 & & 0.0001701751 & *** & MC &   \\
TQMC-wh_n2 vs. ASTRID-ws & abayes & & 50 & 1600 & & 34 & 10 & 6 & & 2.230378e-05 & **** & MC &   \\
TQMC-wh_n2 vs. ASTRID-ws & abayes & & 200 & 200 & & 27 & 7 & 16 & & 0.0003385869 & *** & MC &   \\
TQMC-wh_n2 vs. ASTRID-ws & abayes & & 200 & 400 & & 25 & 13 & 12 & & 0.02814024 & * &  &   \\
TQMC-wh_n2 vs. ASTRID-ws & abayes & & 200 & 800 & & 25 & 13 & 12 & & 0.03383093 & * &  &   \\
TQMC-wh_n2 vs. ASTRID-ws & abayes & & 200 & 1600 & & 29 & 7 & 14 & & 0.00630277 & * &  &   \\
TQMC-wh_n2 vs. ASTRID-ws & abayes & & 500 & 200 & & 24 & 16 & 10 & & 0.1090007 &  &  &   \\
TQMC-wh_n2 vs. ASTRID-ws & abayes & & 500 & 400 & & 23 & 13 & 14 & & 0.1386635 &  &  &   \\
TQMC-wh_n2 vs. ASTRID-ws & abayes & & 500 & 800 & & 22 & 11 & 17 & & 0.02431873 & * &  &   \\
TQMC-wh_n2 vs. ASTRID-ws & abayes & & 500 & 1600 & & 24 & 9 & 17 & & 0.03998109 & * &  &   \\
TQMC-wh_n2 vs. ASTRID-ws & abayes & & 1000 & 200 & & 25 & 9 & 16 & & 0.008141451 & * &  &   \\
TQMC-wh_n2 vs. ASTRID-ws & abayes & & 1000 & 400 & & 21 & 13 & 16 & & 0.07196359 &  &  &   \\
TQMC-wh_n2 vs. ASTRID-ws & abayes & & 1000 & 800 & & 21 & 8 & 21 & & 0.004363813 & ** &  &   \\
TQMC-wh_n2 vs. ASTRID-ws & abayes & & 1000 & 1600 & & 23 & 5 & 22 & & 0.0003702417 & *** & MC &   \\
TQMC-wh_n2 vs. ASTRID-ws & bs & & 50 & 200 & & 30 & 12 & 8 & & 0.004163601 & ** &  &   \\
TQMC-wh_n2 vs. ASTRID-ws & bs & & 50 & 400 & & 31 & 12 & 7 & & 0.002196742 & ** &  &   \\
TQMC-wh_n2 vs. ASTRID-ws & bs & & 50 & 800 & & 31 & 13 & 6 & & 0.002441709 & ** &  &   \\
TQMC-wh_n2 vs. ASTRID-ws & bs & & 50 & 1600 & & 32 & 9 & 9 & & 0.0002453604 & *** & MC &   \\
TQMC-wh_n2 vs. ASTRID-ws & bs & & 200 & 200 & & 20 & 18 & 12 & & 0.5903075 &  &  &   \\
TQMC-wh_n2 vs. ASTRID-ws & bs & & 200 & 400 & & 22 & 14 & 14 & & 0.03486185 & * &  &   \\
TQMC-wh_n2 vs. ASTRID-ws & bs & & 200 & 800 & & 23 & 11 & 16 & & 0.02033664 & * &  &   \\
TQMC-wh_n2 vs. ASTRID-ws & bs & & 200 & 1600 & & 19 & 11 & 20 & & 0.03968453 & * &  &   \\
TQMC-wh_n2 vs. ASTRID-ws & bs & & 500 & 200 & & 17 & 14 & 19 & & 0.6101888 &  &  &   \\
TQMC-wh_n2 vs. ASTRID-ws & bs & & 500 & 400 & & 25 & 9 & 16 & & 0.001592791 & ** &  &   \\
TQMC-wh_n2 vs. ASTRID-ws & bs & & 500 & 800 & & 17 & 12 & 21 & & 0.1690782 &  &  &   \\
TQMC-wh_n2 vs. ASTRID-ws & bs & & 500 & 1600 & & 20 & 9 & 21 & & 0.02918072 & * &  &   \\
TQMC-wh_n2 vs. ASTRID-ws & bs & & 1000 & 200 & & 23 & 18 & 9 & & 1 &  &  &   \\
TQMC-wh_n2 vs. ASTRID-ws & bs & & 1000 & 400 & & 21 & 11 & 18 & & 0.09418824 &  &  &   \\
TQMC-wh_n2 vs. ASTRID-ws & bs & & 1000 & 800 & & 15 & 9 & 26 & & 0.1307372 &  &  &   \\
TQMC-wh_n2 vs. ASTRID-ws & bs & & 1000 & 1600 & & 15 & 12 & 23 & & 0.4097079 &  &  &   \\

TQMC-wh_n2 vs. TQMC-n2-refinedpoly & none & & 50 & 200 & & 38 & 8 & 4 & & 2.767127e-07 & ***** & MC &   \\
TQMC-wh_n2 vs. TQMC-n2-refinedpoly & none & & 50 & 400 & & 31 & 11 & 8 & & 0.0001133465 & *** & MC &   \\
TQMC-wh_n2 vs. TQMC-n2-refinedpoly & none & & 50 & 800 & & 32 & 8 & 10 & & 7.426644e-05 & *** & MC &   \\
TQMC-wh_n2 vs. TQMC-n2-refinedpoly & none & & 50 & 1600 & & 24 & 16 & 10 & & 0.05026649 &  &  &   \\
TQMC-wh_n2 vs. TQMC-n2-refinedpoly & none & & 200 & 200 & & 35 & 11 & 4 & & 1.377439e-06 & ***** & MC &   \\
TQMC-wh_n2 vs. TQMC-n2-refinedpoly & none & & 200 & 400 & & 34 & 8 & 8 & & 6.63259e-06 & **** & MC &   \\
TQMC-wh_n2 vs. TQMC-n2-refinedpoly & none & & 200 & 800 & & 22 & 11 & 17 & & 0.02868212 & * &  &   \\
TQMC-wh_n2 vs. TQMC-n2-refinedpoly & none & & 200 & 1600 & & 29 & 9 & 12 & & 0.0008714751 & ** &  &   \\
TQMC-wh_n2 vs. TQMC-n2-refinedpoly & none & & 500 & 200 & & 29 & 9 & 12 & & 0.0005919595 & ** & MC &   \\
TQMC-wh_n2 vs. TQMC-n2-refinedpoly & none & & 500 & 400 & & 28 & 9 & 13 & & 0.0001411493 & *** & MC &   \\
TQMC-wh_n2 vs. TQMC-n2-refinedpoly & none & & 500 & 800 & & 22 & 9 & 19 & & 0.007677527 & * &  &   \\
TQMC-wh_n2 vs. TQMC-n2-refinedpoly & none & & 500 & 1600 & & 24 & 6 & 20 & & 0.0002278164 & *** & MC &   \\
TQMC-wh_n2 vs. TQMC-n2-refinedpoly & none & & 1000 & 200 & & 26 & 6 & 18 & & 4.012324e-05 & **** & MC &   \\
TQMC-wh_n2 vs. TQMC-n2-refinedpoly & none & & 1000 & 400 & & 29 & 4 & 17 & & 8.018687e-07 & ***** & MC &   \\
TQMC-wh_n2 vs. TQMC-n2-refinedpoly & none & & 1000 & 800 & & 28 & 5 & 17 & & 8.580275e-06 & **** & MC &   \\
TQMC-wh_n2 vs. TQMC-n2-refinedpoly & none & & 1000 & 1600 & & 22 & 5 & 23 & & 0.001707599 & ** &  &   \\

[1] "done"