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

        # Remove ties for wilcox test
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

            # Perform wilcoxon test
            wsr <- wilcox.test(tmp1, 
                               tmp2,
                               paired=TRUE,
                               alternative="two.sided",
                               mu=0, exact=NULL)

            # Report
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
                             supp, "&",
                             as.character(ngen), "&",
                             as.character(nbps), "&",
                             as.character(ntreeqmc), "/",
                             as.character(nother), "/",
                             as.character(ntie), "&",
                             format(wsr$p.value, scientific = TRUE), "&", stars, "&", mc)
                print(writeme)
        } else {
            writeme <- paste("TQMC-wh-n2 vs.", mthd, "&", 
                             supp, "&",
                             as.character(ngen), "&",
                             as.character(nbps), "&",
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


#[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 50 & 200 & 23 / 17 / 10 & 2.951435e-01 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 50 & 400 & 33 / 7 / 10 & 4.793529e-04 & *** & MC"
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 50 & 800 & 27 / 8 / 15 & 1.112766e-04 & *** & MC"
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 50 & 1600 & 23 / 6 / 21 & 9.786145e-04 & ** & "
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 200 & 200 & 30 / 9 / 11 & 5.375914e-03 & * & "
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 200 & 400 & 21 / 14 / 15 & 1.422518e-01 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 200 & 800 & 18 / 11 / 21 & 2.293086e-01 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 200 & 1600 & 22 / 10 / 18 & 2.262861e-02 & * & "
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 500 & 200 & 22 / 15 / 13 & 3.525035e-01 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 500 & 400 & 14 / 22 / 14 & 7.207358e-01 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 500 & 800 & 18 / 11 / 21 & 4.303105e-01 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 500 & 1600 & 21 / 13 / 16 & 1.607812e-01 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 1000 & 200 & 20 / 9 / 21 & 1.403843e-03 & ** & "
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 1000 & 400 & 20 / 6 / 24 & 9.810156e-03 & * & "
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 1000 & 800 & 20 / 7 / 23 & 1.364808e-02 & * & "
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 1000 & 1600 & 19 / 5 / 26 & & NA & NA & NA"

#[1] "TQMC-wh-n2 vs. ASTER-h & bs & 50 & 200 & 22 / 17 / 11 & 3.6032e-01 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & bs & 50 & 400 & 19 / 14 / 17 & 3.014786e-01 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & bs & 50 & 800 & 23 / 18 / 9 & 7.003239e-01 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & bs & 50 & 1600 & 21 / 16 / 13 & 2.202192e-01 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & bs & 200 & 200 & 19 / 17 / 14 & 4.623817e-01 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & bs & 200 & 400 & 21 / 11 / 18 & 7.799515e-02 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & bs & 200 & 800 & 17 / 18 / 15 & 5.335754e-01 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & bs & 200 & 1600 & 17 / 11 / 22 & 2.138936e-01 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & bs & 500 & 200 & 23 / 18 / 9 & 6.237101e-01 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & bs & 500 & 400 & 21 / 14 / 15 & 8.318931e-02 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & bs & 500 & 800 & 17 / 10 / 23 & 3.587861e-01 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & bs & 500 & 1600 & 13 / 8 / 29 & & NA & NA & NA"
#[1] "TQMC-wh-n2 vs. ASTER-h & bs & 1000 & 200 & 16 / 19 / 15 & 5.822796e-01 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & bs & 1000 & 400 & 16 / 10 / 24 & 6.036127e-02 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & bs & 1000 & 800 & 17 / 6 / 27 & & NA & NA & NA"
#[1] "TQMC-wh-n2 vs. ASTER-h & bs & 1000 & 1600 & 13 / 10 / 27 & & NA & NA & NA"

#[1] "TQMC-wh-n2 vs. wASTRID & abayes & 50 & 200 & 30 / 16 / 4 & 2.296252e-02 & * & "
#[1] "TQMC-wh-n2 vs. wASTRID & abayes & 50 & 400 & 33 / 12 / 5 & 1.416331e-03 & ** & "
#[1] "TQMC-wh-n2 vs. wASTRID & abayes & 50 & 800 & 35 / 10 / 5 & 9.071185e-05 & *** & MC"
#[1] "TQMC-wh-n2 vs. wASTRID & abayes & 50 & 1600 & 34 / 10 / 6 & 7.205789e-05 & *** & MC"
#[1] "TQMC-wh-n2 vs. wASTRID & abayes & 200 & 200 & 27 / 7 / 16 & 2.912295e-04 & *** & MC"
#[1] "TQMC-wh-n2 vs. wASTRID & abayes & 200 & 400 & 25 / 13 / 12 & 1.913904e-02 & * & "
#[1] "TQMC-wh-n2 vs. wASTRID & abayes & 200 & 800 & 25 / 13 / 12 & 1.182815e-02 & * & "
#[1] "TQMC-wh-n2 vs. wASTRID & abayes & 200 & 1600 & 29 / 7 / 14 & 6.481961e-03 & * & "
#[1] "TQMC-wh-n2 vs. wASTRID & abayes & 500 & 200 & 24 / 16 / 10 & 3.823234e-02 & * & "
#[1] "TQMC-wh-n2 vs. wASTRID & abayes & 500 & 400 & 23 / 13 / 14 & 1.562814e-01 &  & "
#[1] "TQMC-wh-n2 vs. wASTRID & abayes & 500 & 800 & 22 / 11 / 17 & 3.37693e-02 & * & "
#[1] "TQMC-wh-n2 vs. wASTRID & abayes & 500 & 1600 & 24 / 9 / 17 & 2.002377e-02 & * & "
#[1] "TQMC-wh-n2 vs. wASTRID & abayes & 1000 & 200 & 25 / 9 / 16 & 5.746979e-03 & * & "
#[1] "TQMC-wh-n2 vs. wASTRID & abayes & 1000 & 400 & 21 / 13 / 16 & 3.046085e-02 & * & "
#[1] "TQMC-wh-n2 vs. wASTRID & abayes & 1000 & 800 & 21 / 8 / 21 & 4.159118e-03 & ** & "
#[1] "TQMC-wh-n2 vs. wASTRID & abayes & 1000 & 1600 & 23 / 5 / 22 & 6.004268e-04 & ** & MC"

#[1] "TQMC-wh-n2 vs. wASTRID & bs & 50 & 200 & 30 / 12 / 8 & 3.364415e-03 & ** & "
#[1] "TQMC-wh-n2 vs. wASTRID & bs & 50 & 400 & 31 / 12 / 7 & 3.938716e-03 & ** & "
#[1] "TQMC-wh-n2 vs. wASTRID & bs & 50 & 800 & 31 / 13 / 6 & 1.56284e-03 & ** & "
#[1] "TQMC-wh-n2 vs. wASTRID & bs & 50 & 1600 & 32 / 9 / 9 & 5.733083e-04 & ** & MC"
#[1] "TQMC-wh-n2 vs. wASTRID & bs & 200 & 200 & 20 / 18 / 12 & 3.007714e-01 &  & "
#[1] "TQMC-wh-n2 vs. wASTRID & bs & 200 & 400 & 22 / 14 / 14 & 1.110129e-02 & * & "
#[1] "TQMC-wh-n2 vs. wASTRID & bs & 200 & 800 & 23 / 11 / 16 & 1.713481e-02 & * & "
#[1] "TQMC-wh-n2 vs. wASTRID & bs & 200 & 1600 & 19 / 11 / 20 & 8.739463e-03 & * & "
#[1] "TQMC-wh-n2 vs. wASTRID & bs & 500 & 200 & 17 / 14 / 19 & 4.283361e-01 &  & "
#[1] "TQMC-wh-n2 vs. wASTRID & bs & 500 & 400 & 25 / 9 / 16 & 3.284391e-03 & ** & "
#[1] "TQMC-wh-n2 vs. wASTRID & bs & 500 & 800 & 17 / 12 / 21 & 8.919014e-02 &  & "
#[1] "TQMC-wh-n2 vs. wASTRID & bs & 500 & 1600 & 20 / 9 / 21 & 4.520053e-02 & * & "
#[1] "TQMC-wh-n2 vs. wASTRID & bs & 1000 & 200 & 23 / 18 / 9 & 9.417421e-01 &  & "
#[1] "TQMC-wh-n2 vs. wASTRID & bs & 1000 & 400 & 21 / 11 / 18 & 6.551437e-02 &  & "
#[1] "TQMC-wh-n2 vs. wASTRID & bs & 1000 & 800 & 15 / 9 / 26 & & NA & NA & NA"
#[1] "TQMC-wh-n2 vs. wASTRID & bs & 1000 & 1600 & 15 / 12 / 23 & 6.30483e-01 &  & "
#There were 50 or more warnings (use warnings() to see the first 50)