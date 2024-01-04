#!/usr/bin/env Rscript

# Remove ties
# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0200837

data <- read.csv("../csvs/data-for-testing.csv")
data$NBPS <- as.factor(data$NBPS)
data$NGEN <- as.factor(data$NGEN)
data$SUPP <- as.factor(data$SUPP)
data$REPL <- as.factor(data$REPL)

supps <- c("abayes", "bs")
ngens <- c("50", "200", "500", "1000")
nbpss <- c("200", "400", "800", "1600")
mthds <- c("ASTER-h", "wASTRID")

threshold5 <- 0.05 
threshold4 <- 0.005
threshold3 <- 0.0005
threshold2 <- 0.00005
threshold1 <- 0.000005

ntest <- length(supps) * length(ngens) * length(nbpss) * length(mthds)
threshold_bonferoni <- 0.05 / ntest

for (mthd in mthds) {
    for (supp in supps) {
        for (ngen in ngens) {
            for (nbps in nbpss) {
        
        fndf <- data[data$SUPP == supp, ]
        fndf <- fndf[fndf$NGEN == ngen, ]
        fndf <- fndf[fndf$NBPS == nbps, ]

        # Remove ties
        if (mthd == "ASTER-h") {
            fndf$DIFF <- fndf$ASTERHxSEFN - fndf$TQMCwhn2xSEFN  # positive means TQMC-wh-n2 is better
        } else if (mthd == "wASTRID") {
            fndf$DIFF <- fndf$WASTRIDxSEFN - fndf$TQMCwhn2xSEFN  # positive means TQMC-wh-n2 is better
        } else {
            print("ERROR")
            exit()
        }

        # Count false negatives
        ntie <- sum(fndf$DIFF == 0)
        nother <- sum(fndf$DIFF < 0)
        ntreeqmc <- sum(fndf$DIFF > 0)

        # Remove ties
        fndf$DIFF[fndf$DIFF == 0] <- NA
        xfndf <- fndf[complete.cases(fndf), ]

        if (nrow(xfndf) > 25) {
            tmp1 <- xfndf$TQMCwhn2xSEFN
            if (mthd == "ASTER-h") {
                tmp2 <- xfndf$ASTERHxSEFN
            } else if (mthd == "wASTRID") {
                tmp2 <- xfndf$WASTRIDxSEFN
            } else {
                print("ERROR")
                exit()
            }

            #wsr <- wilcox.test(tmp1, 
            #                   tmp2,
            #                   paired=TRUE,
            #                   alternative="two.sided",
            #                   mu=0, exact=NULL)

            # switch to paired ttest (although similar results)
            wsr <- t.test(tmp1, y = tmp2,
                    alternative="two.sided",
                    mu = 0, paired=TRUE, var.equal=FALSE,
                    conf.level=0.95)

            if (wsr$p.value < 0.000005) {
                stars = "*****"
            } else if (wsr$p.value < 0.00005) {
                stars = "****"
            } else if (wsr$p.value < 0.0005) {
                stars = "***"
            } else if (wsr$p.value < 0.005) {
                stars = "**"
            } else if (wsr$p.value < 0.05) {
                stars = "*"
            } else {
                stars = ""
            }

            if (wsr$p.value < threshold_bonferoni) {
                mc <- "MC"
            } else {
                mc <- ""
            }

            writeme <- paste("TQMC-wh-n2 vs.", mthd, "&",
                             supp, as.character(ngen), as.character(nbps), "&",
                             as.character(ntreeqmc), "/",
                             as.character(nother), "/",
                             as.character(ntie), "&",
                             format(wsr$p.value, scientific = TRUE), "&", stars, "&", mc)
                print(writeme)
        } else {
            writeme <- paste("TQMC-wh-n2 vs.", mthd, "&", 
                             supp, as.character(ngen), as.character(nbps), "&",
                             as.character(ntreeqmc), "/",
                             as.character(nother), "/",
                             as.character(ntie), "&",
                             "& NA & NA & NA")
            print(writeme)
        }
    }

        }
    }
}

#[1] "TQMC-wh-n2 vs. ASTER-h & abayes 50 200 & 23 / 17 / 10 & 2.430175e-01 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes 50 400 & 33 / 7 / 10 & 1.312945e-04 & *** & MC"
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes 50 800 & 27 / 8 / 15 & 1.710143e-05 & **** & MC"
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes 50 1600 & 23 / 6 / 21 & 5.191495e-04 & ** & MC"
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes 200 200 & 30 / 9 / 11 & 1.998263e-03 & ** & "
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes 200 400 & 21 / 14 / 15 & 9.13736e-02 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes 200 800 & 18 / 11 / 21 & 2.558599e-01 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes 200 1600 & 22 / 10 / 18 & 2.07732e-02 & * & "
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes 500 200 & 22 / 15 / 13 & 3.170731e-01 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes 500 400 & 14 / 22 / 14 & 6.58495e-01 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes 500 800 & 18 / 11 / 21 & 2.500752e-01 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes 500 1600 & 21 / 13 / 16 & 1.389488e-01 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes 1000 200 & 20 / 9 / 21 & 2.157428e-03 & ** & "
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes 1000 400 & 20 / 6 / 24 & 1.650162e-02 & * & "
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes 1000 800 & 20 / 7 / 23 & 1.437209e-02 & * & "
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes 1000 1600 & 19 / 5 / 26 & & NA & NA & NA"

#[1] "TQMC-wh-n2 vs. ASTER-h & bs 50 200 & 22 / 17 / 11 & 6.07271e-01 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & bs 50 400 & 19 / 14 / 17 & 1.958497e-01 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & bs 50 800 & 23 / 18 / 9 & 7.546199e-01 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & bs 50 1600 & 21 / 16 / 13 & 2.075663e-01 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & bs 200 200 & 19 / 17 / 14 & 5.951506e-01 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & bs 200 400 & 21 / 11 / 18 & 9.474162e-02 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & bs 200 800 & 17 / 18 / 15 & 4.431427e-01 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & bs 200 1600 & 17 / 11 / 22 & 2.012377e-01 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & bs 500 200 & 23 / 18 / 9 & 9.475904e-01 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & bs 500 400 & 21 / 14 / 15 & 1.008038e-01 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & bs 500 800 & 17 / 10 / 23 & 6.134313e-01 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & bs 500 1600 & 13 / 8 / 29 & & NA & NA & NA"
#[1] "TQMC-wh-n2 vs. ASTER-h & bs 1000 200 & 16 / 19 / 15 & 3.768261e-01 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & bs 1000 400 & 16 / 10 / 24 & 5.715934e-02 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & bs 1000 800 & 17 / 6 / 27 & & NA & NA & NA"
#[1] "TQMC-wh-n2 vs. ASTER-h & bs 1000 1600 & 13 / 10 / 27 & & NA & NA & NA"

#[1] "TQMC-wh-n2 vs. wASTRID & abayes 50 200 & 30 / 16 / 4 & 2.771708e-02 & * & "
#[1] "TQMC-wh-n2 vs. wASTRID & abayes 50 400 & 33 / 12 / 5 & 7.823235e-04 & ** & "
#[1] "TQMC-wh-n2 vs. wASTRID & abayes 50 800 & 35 / 10 / 5 & 6.970537e-05 & *** & MC"
#[1] "TQMC-wh-n2 vs. wASTRID & abayes 50 1600 & 34 / 10 / 6 & 4.27361e-05 & **** & MC"
#[1] "TQMC-wh-n2 vs. wASTRID & abayes 200 200 & 27 / 7 / 16 & 6.07669e-05 & *** & MC"
#[1] "TQMC-wh-n2 vs. wASTRID & abayes 200 400 & 25 / 13 / 12 & 1.152585e-02 & * & "
#[1] "TQMC-wh-n2 vs. wASTRID & abayes 200 800 & 25 / 13 / 12 & 1.539516e-02 & * & "
#[1] "TQMC-wh-n2 vs. wASTRID & abayes 200 1600 & 29 / 7 / 14 & 5.971463e-03 & * & "
#[1] "TQMC-wh-n2 vs. wASTRID & abayes 500 200 & 24 / 16 / 10 & 3.296637e-02 & * & "
#[1] "TQMC-wh-n2 vs. wASTRID & abayes 500 400 & 23 / 13 / 14 & 3.580851e-01 &  & "
#[1] "TQMC-wh-n2 vs. wASTRID & abayes 500 800 & 22 / 11 / 17 & 2.408066e-02 & * & "
#[1] "TQMC-wh-n2 vs. wASTRID & abayes 500 1600 & 24 / 9 / 17 & 1.250169e-02 & * & "
#[1] "TQMC-wh-n2 vs. wASTRID & abayes 1000 200 & 25 / 9 / 16 & 2.44088e-03 & ** & "
#[1] "TQMC-wh-n2 vs. wASTRID & abayes 1000 400 & 21 / 13 / 16 & 2.670654e-02 & * & "
#[1] "TQMC-wh-n2 vs. wASTRID & abayes 1000 800 & 21 / 8 / 21 & 2.860223e-03 & ** & "
#[1] "TQMC-wh-n2 vs. wASTRID & abayes 1000 1600 & 23 / 5 / 22 & 2.101322e-03 & ** & "

#[1] "TQMC-wh-n2 vs. wASTRID & bs 50 200 & 30 / 12 / 8 & 2.88634e-03 & ** & "
#[1] "TQMC-wh-n2 vs. wASTRID & bs 50 400 & 31 / 12 / 7 & 2.159028e-03 & ** & "
#[1] "TQMC-wh-n2 vs. wASTRID & bs 50 800 & 31 / 13 / 6 & 2.954146e-03 & ** & "
#[1] "TQMC-wh-n2 vs. wASTRID & bs 50 1600 & 32 / 9 / 9 & 2.579435e-04 & *** & MC"
#[1] "TQMC-wh-n2 vs. wASTRID & bs 200 200 & 20 / 18 / 12 & 4.381943e-01 &  & "
#[1] "TQMC-wh-n2 vs. wASTRID & bs 200 400 & 22 / 14 / 14 & 7.361387e-03 & * & "
#[1] "TQMC-wh-n2 vs. wASTRID & bs 200 800 & 23 / 11 / 16 & 9.601488e-03 & * & "
#[1] "TQMC-wh-n2 vs. wASTRID & bs 200 1600 & 19 / 11 / 20 & 1.261935e-02 & * & "
#[1] "TQMC-wh-n2 vs. wASTRID & bs 500 200 & 17 / 14 / 19 & 3.839088e-01 &  & "
#[1] "TQMC-wh-n2 vs. wASTRID & bs 500 400 & 25 / 9 / 16 & 4.802092e-03 & ** & "
#[1] "TQMC-wh-n2 vs. wASTRID & bs 500 800 & 17 / 12 / 21 & 1.474068e-01 &  & "
#[1] "TQMC-wh-n2 vs. wASTRID & bs 500 1600 & 20 / 9 / 21 & 3.530596e-02 & * & "
#[1] "TQMC-wh-n2 vs. wASTRID & bs 1000 200 & 23 / 18 / 9 & 8.643007e-01 &  & "
#[1] "TQMC-wh-n2 vs. wASTRID & bs 1000 400 & 21 / 11 / 18 & 6.22272e-02 &  & "
#[1] "TQMC-wh-n2 vs. wASTRID & bs 1000 800 & 15 / 9 / 26 & & NA & NA & NA"
#[1] "TQMC-wh-n2 vs. wASTRID & bs 1000 1600 & 15 / 12 / 23 & 5.564763e-01 &  & "
