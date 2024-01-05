#!/usr/bin/env Rscript

# MANY NORMALITY FAILURES SO USE WILCOXON INSTEAD!!!

data <- read.csv("../csvs/data-for-testing.csv")
data$NBPS <- as.factor(data$NBPS)
data$NGEN <- as.factor(data$NGEN)
data$SUPP <- as.factor(data$SUPP)
data$REPL <- as.factor(data$REPL)

supps <- c("abayes", "bs")
ngens <- c("50", "200", "500", "1000")
nbpss <- c("200", "400", "800", "1600")
mthds <- c("ASTER-h", "wASTRID", "")

threshold0 <- 0.01
ntests <- length(supps) * length(ngens) * length(nbpss) * length(mthds)
threshold_bonferoni <- threshold0 / ntests

writeme <- paste("alpha =", as.character(threshold0))
print(writeme)
print(paste("bonferroni alpha =", as.character(threshold_bonferoni), 
            "=", as.character(threshold0), "/", as.character(ntests)))

for (mthd in mthds) {
    for (supp in supps) {
        for (ngen in ngens) {
            for (nbps in nbpss) {
        
        df <- data[data$SUPP == supp, ]
        xdf <- df[df$NGEN == ngen, ]
        xdf <- df[df$NBPS == nbps, ]

        # Find ties
        mthd1 <- xdf$TQMCwhn2xSEF
        if (mthd == "ASTER-h") {
            mthd2 <- xdf$ASTERHxSEFN
        } else if (mthd == "wASTRID") {
            mthd2 <- xdf$WASTRIDxSEFN
        } else {
            print("ERROR")
            exit()
        }

        # Count ties
        diff <- mthd2 - mthd1
        ntie <- sum(diff == 0)
        nother <- sum(diff < 0)
        ntreeqmc <- sum(diff > 0)
    
        tt <- t.test(mthd1, y = mthd2,
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
            # To be safe, check for normality
            #ggqqplot(diff)
            #hist(diff)
            st <- shapiro.test(diff)

            stars = ""
            if (tt$p.value < 0.000001) {
                stars = "*****"
            } else if (tt$p.value < 0.00001) {
                stars = "****"
            } else if (tt$p.value < 0.0001) {
                stars = "***"
            } else if (tt$p.value < 0.001) {
                stars = "**"
            } else if (tt$p.value < threshold0) {
                stars = "*"
            }
            mc <- ""
            if (tt$p.value < threshold_bonferoni) {
                mc <- "MC"
            }
            note <- ""
            if ((mean(diff) < 0) & (tt$p.value < threshold0)) {
                note <- paste("(", mthd, "better)")
            }
            if (st$p.value < 0.05) {
                note <- paste(note, "<- Failed normality test!")
            }
        }

        writeme <- paste("TQMC-wh-n2 vs.", mthd, "&",
                         supp, "&",
                         as.character(ngen), "&",
                         as.character(nbps), "&",
                         as.character(ntreeqmc), "/",
                         as.character(nother), "/",
                         as.character(ntie), "&",
                         format(tt$p.value, scientific = TRUE),
                         "&", stars, "&", mc,
                         note)
        print(writeme)

            }
        }
    }
}
