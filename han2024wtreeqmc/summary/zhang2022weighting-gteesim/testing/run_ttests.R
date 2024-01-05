#!/usr/bin/env Rscript

data <- read.csv("../csvs/data-for-testing.csv")
data$NBPS <- as.factor(data$NBPS)
data$NGEN <- as.factor(data$NGEN)
data$SUPP <- as.factor(data$SUPP)
data$REPL <- as.factor(data$REPL)

supps <- c("abayes", "bs")
ngens <- c("50", "200", "500", "1000")
nbpss <- c("200", "400", "800", "1600")
mthds <- c("ASTER-h", "wASTRID")

ntest <- length(supps) * length(ngens) * length(nbpss) * length(mthds)
threshold_bonferoni <- 0.05 / ntest

for (mthd in mthds) {
    for (supp in supps) {
        for (ngen in ngens) {
            for (nbps in nbpss) {
        
        fndf <- data[data$SUPP == supp, ]
        fndf <- fndf[fndf$NGEN == ngen, ]
        fndf <- fndf[fndf$NBPS == nbps, ]

        # Find ties
        if (mthd == "ASTER-h") {
            fndf$DIFF <- fndf$ASTERHxSEFN - fndf$TQMCwhn2xSEFN  # positive means TQMC-wh-n2 is better
        } else if (mthd == "wASTRID") {
            fndf$DIFF <- fndf$WASTRIDxSEFN - fndf$TQMCwhn2xSEFN  # positive means TQMC-wh-n2 is better
        } else {
            print("ERROR")
            exit()
        }

        # Count ties
        ntie <- sum(fndf$DIFF == 0)
        nother <- sum(fndf$DIFF < 0)
        ntreeqmc <- sum(fndf$DIFF > 0)

        # Run paired t-test
        #tmp1 <- fndf$TQMCwhn2xSEFN
        tmp1 <- fndf$TQMCwhn2xSERF
        if (mthd == "ASTER-h") {
            #tmp2 <- fndf$ASTERHxSEFN
            tmp2 <- fndf$ASTERHxSERF
        } else if (mthd == "wASTRID") {
            #tmp2 <- fndf$WASTRIDxSEFN
            tmp2 <- fndf$WASTRIDxSERF
        } else {
            print("ERROR")
            exit()
        }
    
        tt <- t.test(tmp1, y = tmp2,
                     alternative="two.sided",
                     mu = 0, paired=TRUE, var.equal=FALSE,
                     conf.level=0.95)

        # Report
        if (ntreeqmc + nother == 0) {
            # all ties!
            tt$p.value <- "NA"
            stars = ""
            mc <- ""
        } else {
            if (tt$p.value < 0.000005) {
                stars = "*****"
            } else if (tt$p.value < 0.00005) {
                stars = "****"
            } else if (tt$p.value < 0.0005) {
                stars = "***"
            } else if (tt$p.value < 0.005) {
                stars = "**"
            } else if (tt$p.value < 0.05) {
                stars = "*"
            } else {
                stars = ""
            }
            if (tt$p.value < threshold_bonferoni) {
                mc <- "MC"
            } else {
                mc <- ""
            }
        }

        if (tt$p.value < threshold_bonferoni) {
            mc <- "MC"
        } else {
            mc <- ""
        }

        writeme <- paste("TQMC-wh-n2 vs.", mthd, "&",
                         supp, "&",
                         as.character(ngen), "&",
                         as.character(nbps), "&",
                         as.character(ntreeqmc), "/",
                         as.character(nother), "/",
                         as.character(ntie), "&",
                         format(tt$p.value, scientific = TRUE), "&", stars, "&", mc)
        print(writeme)

            }
        }
    }
}

#[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 50 & 200 & 23 / 17 / 10 & 2.420411e-01 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 50 & 400 & 33 / 7 / 10 & 1.703796e-04 & *** & MC"
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 50 & 800 & 27 / 8 / 15 & 3.588418e-05 & **** & MC"
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 50 & 1600 & 23 / 6 / 21 & 8.335666e-04 & ** & "
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 200 & 200 & 30 / 9 / 11 & 2.22843e-03 & ** & "
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 200 & 400 & 21 / 14 / 15 & 9.133392e-02 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 200 & 800 & 18 / 11 / 21 & 2.528768e-01 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 200 & 1600 & 22 / 10 / 18 & 2.175882e-02 & * & "
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 500 & 200 & 22 / 15 / 13 & 3.153366e-01 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 500 & 400 & 14 / 22 / 14 & 6.566915e-01 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 500 & 800 & 18 / 11 / 21 & 2.471621e-01 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 500 & 1600 & 21 / 13 / 16 & 1.382333e-01 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 1000 & 200 & 20 / 9 / 21 & 2.81871e-03 & ** & "
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 1000 & 400 & 20 / 6 / 24 & 1.814899e-02 & * & "
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 1000 & 800 & 20 / 7 / 23 & 1.586451e-02 & * & "
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 1000 & 1600 & 19 / 5 / 26 & 7.19011e-03 & * & "

#[1] "TQMC-wh-n2 vs. ASTER-h & bs & 50 & 200 & 22 / 17 / 11 & 6.058411e-01 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & bs & 50 & 400 & 19 / 14 / 17 & 1.94325e-01 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & bs & 50 & 800 & 23 / 18 / 9 & 7.538395e-01 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & bs & 50 & 1600 & 21 / 16 / 13 & 2.064352e-01 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & bs & 200 & 200 & 19 / 17 / 14 & 5.931409e-01 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & bs & 200 & 400 & 21 / 11 / 18 & 9.463517e-02 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & bs & 200 & 800 & 17 / 18 / 15 & 4.407162e-01 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & bs & 200 & 1600 & 17 / 11 / 22 & 1.988249e-01 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & bs & 500 & 200 & 23 / 18 / 9 & 9.474128e-01 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & bs & 500 & 400 & 21 / 14 / 15 & 1.006376e-01 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & bs & 500 & 800 & 17 / 10 / 23 & 6.091112e-01 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & bs & 500 & 1600 & 13 / 8 / 29 & 2.632516e-01 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & bs & 1000 & 200 & 16 / 19 / 15 & 3.744971e-01 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & bs & 1000 & 400 & 16 / 10 / 24 & 5.809428e-02 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & bs & 1000 & 800 & 17 / 6 / 27 & 1.370959e-01 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & bs & 1000 & 1600 & 13 / 10 / 27 & 5.969924e-01 &  & "

#[1] "TQMC-wh-n2 vs. wASTRID & abayes & 50 & 200 & 30 / 16 / 4 & 2.78601e-02 & * & "
#[1] "TQMC-wh-n2 vs. wASTRID & abayes & 50 & 400 & 33 / 12 / 5 & 8.341777e-04 & ** & "
#[1] "TQMC-wh-n2 vs. wASTRID & abayes & 50 & 800 & 35 / 10 / 5 & 8.017552e-05 & *** & MC"
#[1] "TQMC-wh-n2 vs. wASTRID & abayes & 50 & 1600 & 34 / 10 / 6 & 5.193572e-05 & *** & MC"
#[1] "TQMC-wh-n2 vs. wASTRID & abayes & 200 & 200 & 27 / 7 / 16 & 1.102826e-04 & *** & MC"
#[1] "TQMC-wh-n2 vs. wASTRID & abayes & 200 & 400 & 25 / 13 / 12 & 1.2038e-02 & * & "
#[1] "TQMC-wh-n2 vs. wASTRID & abayes & 200 & 800 & 25 / 13 / 12 & 1.593242e-02 & * & "
#[1] "TQMC-wh-n2 vs. wASTRID & abayes & 200 & 1600 & 29 / 7 / 14 & 6.494759e-03 & * & "
#[1] "TQMC-wh-n2 vs. wASTRID & abayes & 500 & 200 & 24 / 16 / 10 & 3.336011e-02 & * & "
#[1] "TQMC-wh-n2 vs. wASTRID & abayes & 500 & 400 & 23 / 13 / 14 & 3.560201e-01 &  & "
#[1] "TQMC-wh-n2 vs. wASTRID & abayes & 500 & 800 & 22 / 11 / 17 & 2.497082e-02 & * & "
#[1] "TQMC-wh-n2 vs. wASTRID & abayes & 500 & 1600 & 24 / 9 / 17 & 1.336366e-02 & * & "
#[1] "TQMC-wh-n2 vs. wASTRID & abayes & 1000 & 200 & 25 / 9 / 16 & 2.878219e-03 & ** & "
#[1] "TQMC-wh-n2 vs. wASTRID & abayes & 1000 & 400 & 21 / 13 / 16 & 2.750445e-02 & * & "
#[1] "TQMC-wh-n2 vs. wASTRID & abayes & 1000 & 800 & 21 / 8 / 21 & 3.612783e-03 & ** & "
#[1] "TQMC-wh-n2 vs. wASTRID & abayes & 1000 & 1600 & 23 / 5 / 22 & 2.817342e-03 & ** & "

#[1] "TQMC-wh-n2 vs. wASTRID & bs & 50 & 200 & 30 / 12 / 8 & 3.070153e-03 & ** & "
#[1] "TQMC-wh-n2 vs. wASTRID & bs & 50 & 400 & 31 / 12 / 7 & 2.294484e-03 & ** & "
#[1] "TQMC-wh-n2 vs. wASTRID & bs & 50 & 800 & 31 / 13 / 6 & 3.086169e-03 & ** & "
#[1] "TQMC-wh-n2 vs. wASTRID & bs & 50 & 1600 & 32 / 9 / 9 & 3.10809e-04 & *** & MC"
#[1] "TQMC-wh-n2 vs. wASTRID & bs & 200 & 200 & 20 / 18 / 12 & 4.364106e-01 &  & "
#[1] "TQMC-wh-n2 vs. wASTRID & bs & 200 & 400 & 22 / 14 / 14 & 7.923524e-03 & * & "
#[1] "TQMC-wh-n2 vs. wASTRID & bs & 200 & 800 & 23 / 11 / 16 & 1.034169e-02 & * & "
#[1] "TQMC-wh-n2 vs. wASTRID & bs & 200 & 1600 & 19 / 11 / 20 & 1.374946e-02 & * & "
#[1] "TQMC-wh-n2 vs. wASTRID & bs & 500 & 200 & 17 / 14 / 19 & 3.805485e-01 &  & "
#[1] "TQMC-wh-n2 vs. wASTRID & bs & 500 & 400 & 25 / 9 / 16 & 5.390519e-03 & * & "
#[1] "TQMC-wh-n2 vs. wASTRID & bs & 500 & 800 & 17 / 12 / 21 & 1.461453e-01 &  & "
#[1] "TQMC-wh-n2 vs. wASTRID & bs & 500 & 1600 & 20 / 9 / 21 & 3.644464e-02 & * & "
#[1] "TQMC-wh-n2 vs. wASTRID & bs & 1000 & 200 & 23 / 18 / 9 & 8.638486e-01 &  & "
#[1] "TQMC-wh-n2 vs. wASTRID & bs & 1000 & 400 & 21 / 11 / 18 & 6.269462e-02 &  & "
#[1] "TQMC-wh-n2 vs. wASTRID & bs & 1000 & 800 & 15 / 9 / 26 & 4.979648e-02 & * & "
#[1] "TQMC-wh-n2 vs. wASTRID & bs & 1000 & 1600 & 15 / 12 / 23 & 5.518366e-01 &  & "
