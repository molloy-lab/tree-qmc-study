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
