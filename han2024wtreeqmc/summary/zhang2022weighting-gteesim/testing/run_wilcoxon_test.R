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
mthds <- c("ASTER-h", "wASTRID", "TREEQMC-n2-orig")  # with polys refined

threshold0 <- 0.05
ntests <- length(ngens) * length(nbpss) * 5
threshold_bonferoni <- threshold0 / ntests

writeme <- paste("alpha =", as.character(threshold0))
print(writeme)
print(paste("bonferroni alpha =", as.character(threshold_bonferoni), 
            "=", as.character(threshold0), "/", as.character(ntests)))

for (mthd in mthds) {

    if (mthd == "TREEQMC-n2-orig") {
        supps <- c("abayes")
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
        if (mthd == "ASTER-h") {
            mthd2 <- df$ASTERHxSERF #FN
        } else if (mthd == "wASTRID") {
            mthd2 <- df$WASTRIDxSERF #FN
        } else if (mthd == "TREEQMC-n2-orig") {
            mthd2 <- df$TQMCwnn2xSERF #FN
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

        writeme <- paste("TQMC-wh-n2 vs.", mthd, "&",
                         supp, "&",
                         as.character(ngen), "&",
                             as.character(nbps), "&",
                             as.character(ntreeqmc), "/",
                             as.character(nother), "/",
                             as.character(ntie), "&",
                             format(xwsr, scientific = TRUE), "&", 
                             format(ywsr, scientific = TRUE), "&", 
                             stars, "&", mc)
        print(writeme)

            }
        }
    }
}

print("done")
exit()

[1] "alpha = 0.05"
[1] "bonferroni alpha = 0.000625 = 0.05 / 80"
[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 50 & 200 & 23 / 17 / 10 & 0.3341381 & 0.3238586 &  & "
[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 50 & 400 & 33 / 7 / 10 & 0.000387488 & 0.0001084018 & *** & MC"
[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 50 & 800 & 27 / 8 / 15 & 0.0001895597 & 0.0002738126 & *** & MC"
[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 50 & 1600 & 23 / 6 / 21 & 0.0002571382 & 0.0004458502 & *** & MC"
[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 200 & 200 & 30 / 9 / 11 & 0.006861866 & 0.002451796 & * & "
[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 200 & 400 & 21 / 14 / 15 & 0.1247336 & 0.1534584 &  & "
[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 200 & 800 & 18 / 11 / 21 & 0.139051 & 0.1550322 &  & "
[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 200 & 1600 & 22 / 10 / 18 & 0.02106874 & 0.02215265 & * & "
[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 500 & 200 & 22 / 15 / 13 & 0.7603878 & 0.5231943 &  & "
[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 500 & 400 & 14 / 22 / 14 & 0.529488 & 0.3493314 &  & "
[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 500 & 800 & 18 / 11 / 21 & 0.5209156 & 0.30146 &  & "
[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 500 & 1600 & 21 / 13 / 16 & 0.1982383 & 0.1719252 &  & "
[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 1000 & 200 & 20 / 9 / 21 & 0.0004676618 & 0.006101582 & * & "
[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 1000 & 400 & 20 / 6 / 24 & 0.003161758 & 0.003370821 & ** & "
[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 1000 & 800 & 20 / 7 / 23 & 0.007444784 & 0.007649198 & * & "
[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 1000 & 1600 & 19 / 5 / 26 & 0.0009840727 & 0.001789212 & ** & "
[1] "TQMC-wh-n2 vs. ASTER-h & bs & 50 & 200 & 22 / 17 / 11 & 0.3399764 & 0.3525142 &  & "
[1] "TQMC-wh-n2 vs. ASTER-h & bs & 50 & 400 & 19 / 14 / 17 & 0.4983073 & 0.4291406 &  & "
[1] "TQMC-wh-n2 vs. ASTER-h & bs & 50 & 800 & 23 / 18 / 9 & 0.8386277 & 0.7120619 &  & "
[1] "TQMC-wh-n2 vs. ASTER-h & bs & 50 & 1600 & 21 / 16 / 13 & 0.1315586 & 0.1983822 &  & "
[1] "TQMC-wh-n2 vs. ASTER-h & bs & 200 & 200 & 19 / 17 / 14 & 0.5360853 & 0.6023837 &  & "
[1] "TQMC-wh-n2 vs. ASTER-h & bs & 200 & 400 & 21 / 11 / 18 & 0.06322916 & 0.0613162 &  & "
[1] "TQMC-wh-n2 vs. ASTER-h & bs & 200 & 800 & 17 / 18 / 15 & 0.3616252 & 0.6384818 &  & "
[1] "TQMC-wh-n2 vs. ASTER-h & bs & 200 & 1600 & 17 / 11 / 22 & 0.2169209 & 0.2262515 &  & "
[1] "TQMC-wh-n2 vs. ASTER-h & bs & 500 & 200 & 23 / 18 / 9 & 0.7983978 & 0.6835137 &  & "
[1] "TQMC-wh-n2 vs. ASTER-h & bs & 500 & 400 & 21 / 14 / 15 & 0.07098527 & 0.112634 &  & "
[1] "TQMC-wh-n2 vs. ASTER-h & bs & 500 & 800 & 17 / 10 / 23 & 0.3459887 & 0.2246011 &  & "
[1] "TQMC-wh-n2 vs. ASTER-h & bs & 500 & 1600 & 13 / 8 / 29 & 0.1746035 & 0.2260199 &  & "
[1] "TQMC-wh-n2 vs. ASTER-h & bs & 1000 & 200 & 16 / 19 / 15 & 0.6453902 & 0.6213272 &  & "
[1] "TQMC-wh-n2 vs. ASTER-h & bs & 1000 & 400 & 16 / 10 / 24 & 0.03263813 & 0.1127932 &  & "
[1] "TQMC-wh-n2 vs. ASTER-h & bs & 1000 & 800 & 17 / 6 / 27 & 0.08854365 & 0.03046513 &  & "
[1] "TQMC-wh-n2 vs. ASTER-h & bs & 1000 & 1600 & 13 / 10 / 27 & 0.2587008 & 0.4147711 &  & "

[1] "TQMC-wh-n2 vs. wASTRID & abayes & 50 & 200 & 30 / 16 / 4 & 0.02684553 & 0.02592068 & * & "
[1] "TQMC-wh-n2 vs. wASTRID & abayes & 50 & 400 & 33 / 12 / 5 & 0.001025001 & 0.0008868985 & ** & "
[1] "TQMC-wh-n2 vs. wASTRID & abayes & 50 & 800 & 35 / 10 / 5 & 0.0001701751 & 0.0001248617 & *** & MC"
[1] "TQMC-wh-n2 vs. wASTRID & abayes & 50 & 1600 & 34 / 10 / 6 & 1.910546e-05 & 2.230378e-05 & **** & MC"
[1] "TQMC-wh-n2 vs. wASTRID & abayes & 200 & 200 & 27 / 7 / 16 & 0.0003385869 & 0.0002484249 & *** & MC"
[1] "TQMC-wh-n2 vs. wASTRID & abayes & 200 & 400 & 25 / 13 / 12 & 0.02555467 & 0.02814024 & * & "
[1] "TQMC-wh-n2 vs. wASTRID & abayes & 200 & 800 & 25 / 13 / 12 & 0.03361491 & 0.03383093 & * & "
[1] "TQMC-wh-n2 vs. wASTRID & abayes & 200 & 1600 & 29 / 7 / 14 & 0.00630277 & 0.001176633 & * & "
[1] "TQMC-wh-n2 vs. wASTRID & abayes & 500 & 200 & 24 / 16 / 10 & 0.09150247 & 0.1090007 &  & "
[1] "TQMC-wh-n2 vs. wASTRID & abayes & 500 & 400 & 23 / 13 / 14 & 0.1386635 & 0.1085876 &  & "
[1] "TQMC-wh-n2 vs. wASTRID & abayes & 500 & 800 & 22 / 11 / 17 & 0.01520404 & 0.02431873 & * & "
[1] "TQMC-wh-n2 vs. wASTRID & abayes & 500 & 1600 & 24 / 9 / 17 & 0.03998109 & 0.01614823 & * & "
[1] "TQMC-wh-n2 vs. wASTRID & abayes & 1000 & 200 & 25 / 9 / 16 & 0.008141451 & 0.005204419 & * & "
[1] "TQMC-wh-n2 vs. wASTRID & abayes & 1000 & 400 & 21 / 13 / 16 & 0.04045433 & 0.07196359 &  & "
[1] "TQMC-wh-n2 vs. wASTRID & abayes & 1000 & 800 & 21 / 8 / 21 & 0.001367025 & 0.004363813 & ** & "
[1] "TQMC-wh-n2 vs. wASTRID & abayes & 1000 & 1600 & 23 / 5 / 22 & 0.0003702417 & 0.0002854466 & *** & MC"
[1] "TQMC-wh-n2 vs. wASTRID & bs & 50 & 200 & 30 / 12 / 8 & 0.004163601 & 0.003493408 & ** & "
[1] "TQMC-wh-n2 vs. wASTRID & bs & 50 & 400 & 31 / 12 / 7 & 0.002196742 & 0.001926339 & ** & "
[1] "TQMC-wh-n2 vs. wASTRID & bs & 50 & 800 & 31 / 13 / 6 & 0.002441709 & 0.002425262 & ** & "
[1] "TQMC-wh-n2 vs. wASTRID & bs & 50 & 1600 & 32 / 9 / 9 & 0.0002453604 & 0.0001671805 & *** & MC"
[1] "TQMC-wh-n2 vs. wASTRID & bs & 200 & 200 & 20 / 18 / 12 & 0.5289875 & 0.5903075 &  & "
[1] "TQMC-wh-n2 vs. wASTRID & bs & 200 & 400 & 22 / 14 / 14 & 0.01087662 & 0.03486185 & * & "
[1] "TQMC-wh-n2 vs. wASTRID & bs & 200 & 800 & 23 / 11 / 16 & 0.01631613 & 0.02033664 & * & "
[1] "TQMC-wh-n2 vs. wASTRID & bs & 200 & 1600 & 19 / 11 / 20 & 0.008662751 & 0.03968453 & * & "
[1] "TQMC-wh-n2 vs. wASTRID & bs & 500 & 200 & 17 / 14 / 19 & 0.6101888 & 0.589102 &  & "
[1] "TQMC-wh-n2 vs. wASTRID & bs & 500 & 400 & 25 / 9 / 16 & 0.0009101708 & 0.001592791 & ** & "
[1] "TQMC-wh-n2 vs. wASTRID & bs & 500 & 800 & 17 / 12 / 21 & 0.0654503 & 0.1690782 &  & "
[1] "TQMC-wh-n2 vs. wASTRID & bs & 500 & 1600 & 20 / 9 / 21 & 0.02840782 & 0.02918072 & * & "
[1] "TQMC-wh-n2 vs. wASTRID & bs & 1000 & 200 & 23 / 18 / 9 & 1 & 0.8306194 &  & "
[1] "TQMC-wh-n2 vs. wASTRID & bs & 1000 & 400 & 21 / 11 / 18 & 0.09418824 & 0.07568016 &  & "
[1] "TQMC-wh-n2 vs. wASTRID & bs & 1000 & 800 & 15 / 9 / 26 & 0.05144358 & 0.1307372 &  & "
[1] "TQMC-wh-n2 vs. wASTRID & bs & 1000 & 1600 & 15 / 12 / 23 & 0.2560448 & 0.4097079 &  & "

[1] "TQMC-wh-n2 vs. TREEQMC-n2-orig & abayes & 50 & 200 & 38 / 8 / 4 & 2.657253e-07 & 2.767127e-07 & ***** & MC"
[1] "TQMC-wh-n2 vs. TREEQMC-n2-orig & abayes & 50 & 400 & 31 / 11 / 8 & 6.426154e-05 & 0.0001133465 & *** & MC"
[1] "TQMC-wh-n2 vs. TREEQMC-n2-orig & abayes & 50 & 800 & 32 / 8 / 10 & 7.426644e-05 & 5.095752e-05 & *** & MC"
[1] "TQMC-wh-n2 vs. TREEQMC-n2-orig & abayes & 50 & 1600 & 24 / 16 / 10 & 0.02891075 & 0.05026649 &  & "
[1] "TQMC-wh-n2 vs. TREEQMC-n2-orig & abayes & 200 & 200 & 35 / 11 / 4 & 7.055397e-07 & 1.377439e-06 & ***** & MC"
[1] "TQMC-wh-n2 vs. TREEQMC-n2-orig & abayes & 200 & 400 & 34 / 8 / 8 & 6.63259e-06 & 6.460657e-06 & **** & MC"
[1] "TQMC-wh-n2 vs. TREEQMC-n2-orig & abayes & 200 & 800 & 22 / 11 / 17 & 0.02142979 & 0.02868212 & * & "
[1] "TQMC-wh-n2 vs. TREEQMC-n2-orig & abayes & 200 & 1600 & 29 / 9 / 12 & 0.0008714751 & 0.0006293832 & ** & "
[1] "TQMC-wh-n2 vs. TREEQMC-n2-orig & abayes & 500 & 200 & 29 / 9 / 12 & 0.0005919595 & 0.0004866494 & ** & MC"
[1] "TQMC-wh-n2 vs. TREEQMC-n2-orig & abayes & 500 & 400 & 28 / 9 / 13 & 5.13956e-05 & 0.0001411493 & *** & MC"
[1] "TQMC-wh-n2 vs. TREEQMC-n2-orig & abayes & 500 & 800 & 22 / 9 / 19 & 0.004579221 & 0.007677527 & * & "
[1] "TQMC-wh-n2 vs. TREEQMC-n2-orig & abayes & 500 & 1600 & 24 / 6 / 20 & 0.0001091119 & 0.0002278164 & *** & MC"
[1] "TQMC-wh-n2 vs. TREEQMC-n2-orig & abayes & 1000 & 200 & 26 / 6 / 18 & 1.190277e-05 & 4.012324e-05 & **** & MC"
[1] "TQMC-wh-n2 vs. TREEQMC-n2-orig & abayes & 1000 & 400 & 29 / 4 / 17 & 3.792811e-07 & 8.018687e-07 & ***** & MC"
[1] "TQMC-wh-n2 vs. TREEQMC-n2-orig & abayes & 1000 & 800 & 28 / 5 / 17 & 7.05081e-06 & 8.580275e-06 & **** & MC"
[1] "TQMC-wh-n2 vs. TREEQMC-n2-orig & abayes & 1000 & 1600 & 22 / 5 / 23 & 0.001707599 & 0.0007659644 & ** & "
[1] "done"