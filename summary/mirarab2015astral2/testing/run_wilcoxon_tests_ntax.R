#!/usr/bin/env Rscript

# Remove ties
# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0200837

data <- read.csv("../csvs/data-for-testing.csv")
data$NTAX <- as.factor(data$NTAX)
data$STRHT <- as.factor(data$STRHT)
data$SRATE <- as.factor(data$SRATE)
data$GTRE <- as.factor(data$GTRE)
data$NGEN <- as.factor(data$NGEN)
data$REPL <- as.factor(data$REPL)

data <- subset(data, STRHT == "2000000")
data <- subset(data, SRATE == "1e-06")
data <- subset(data, GTRE == "estimatedgenetre")
data <- subset(data, NGEN == 1000)

ntaxs <- c(10, 50, 100, 200, 500, 1000)
mthds <- c("FASTRAL", "ASTRAL3")

total <- length(ntaxs) * length(mthds)

threshold2 <- 0.05
threshold1 <- 0.005
threshold0 <- threshold2 / total

writeme <- paste("* indicates p <",  as.character(threshold2))
print(writeme)
writeme <- paste("** indicates p <",  as.character(threshold1))
print(writeme)
writeme <- paste("Note: p <",  as.character(threshold0), "if correcting for", as.character(total), "comparisons")
print(writeme)


for (ntax in ntaxs) {
    ydf <- data[data$NTAX == ntax, ]

    for (mthd in mthds) {

        # Remove ties
        if (mthd == "FASTRAL") {
            ydf$DIFF <- ydf$FASTRALxSEFN - ydf$TQMCn2xSEFN  # positive means TQMCn2 is better
        } else {
            ydf$DIFF <- ydf$ASTRAL3xSEFN - ydf$TQMCn2xSEFN  # positive means TQMCn2 is better
        }

		#if (0 %in% ydf$DIFF) {
		#	print("Found a tie!")
		#}

        ntreeqmc <- sum(ydf$DIFF > 0)
		nother <- sum(ydf$DIFF < 0)
        ntie <- sum(ydf$DIFF == 0)
		
        avgdiffwtie <- mean(ydf$DIFF)
        vardiffwtie <- sd(ydf$DIFF)
        meddiffwtie <- median(ydf$DIFF)

        # Remove ties
		ydf$DIFF[ydf$DIFF == 0] <- NA
 		zdf <- ydf[complete.cases(ydf), ]

        avgdiff <- mean(zdf$DIFF)
        meddiff <- median(zdf$DIFF)

        if (nrow(zdf) > 25) {
            tmp1 <- zdf$TQMCn2xSEFN
            if (mthd == "FASTRAL") {
                tmp2 <- zdf$FASTRALxSEFN
            } else {
                tmp2 <- zdf$ASTRAL3xSEFN
            }

            wsr <- wilcox.test(tmp1, 
                               tmp2,
                               paired=TRUE,
                               alternative="two.sided",
                               mu=0, exact=NULL)
            #print(wsr)

            if (wsr$p.value < threshold1) {
                writeme <- paste(as.character(ntax), "&",
                                 "TQMC-n2 vs.", mthd, "&", 
                                 as.character(ntreeqmc), "/",
                                 as.character(nother), "/",
                                 as.character(ntie), "&",
                                 as.character(avgdiffwtie), "$ pm $",
                                 as.character(vardiffwtie), "&",
                                 format(wsr$p.value, scientific = TRUE), "& ** \\")
                print(writeme)
            } else if (wsr$p.value < threshold2) {
                writeme <- paste(as.character(ntax), "&",
                                 "TQMC-n2 vs.", mthd, "&", 
                                 as.character(ntreeqmc), "/",
                                 as.character(nother), "/",
                                 as.character(ntie), "&",
                                 as.character(avgdiffwtie), "$ pm $",
                                 as.character(vardiffwtie), "&",
                                 format(wsr$p.value, scientific = TRUE), "& * \\")
                print(writeme)
            } else {
                writeme <- paste(as.character(ntax), "&",
                                 "TQMC-n2 vs.", mthd, "&", 
                                 as.character(ntreeqmc), "/",
                                 as.character(nother), "/",
                                 as.character(ntie), "&",
                                 as.character(avgdiffwtie), "$ pm $",
                                 as.character(vardiffwtie), "&",
                                 format(wsr$p.value, scientific = FALSE), "\\")
                print(writeme)
            }
        } else {
            writeme <- paste(as.character(ntax), "&",
                             "TQMC-n2 vs.", mthd, "&", 
                             as.character(ntreeqmc), "/",
                             as.character(nother), "/",
                             as.character(ntie), "&",
                             as.character(avgdiffwtie), "$ pm $",
                             as.character(vardiffwtie), "& NA & NA \\")
            print(writeme)
        }
    }
}

#[1] "* indicates p < 0.05"
#[1] "** indicates p < 0.005"
#[1] "Note: p < 0.00416666666666667 if correcting for 12 comparisons"
#
#[1] "10 & TQMC-n2 vs. FASTRAL & 0 / 0 / 49 & 0 $ pm $ 0 & NA & NA \\"
#[1] "10 & TQMC-n2 vs. ASTRAL3 & 0 / 0 / 49 & 0 $ pm $ 0 & NA & NA \\"
#
#[1] "50 & TQMC-n2 vs. FASTRAL & 5 / 3 / 40 & 0.0416666666666667 $ pm $ 0.410414079086058 & NA & NA \\"
#[1] "50 & TQMC-n2 vs. ASTRAL3 & 5 / 3 / 40 & 0.0416666666666667 $ pm $ 0.410414079086058 & NA & NA \\"
#
#[1] "100 & TQMC-n2 vs. FASTRAL & 9 / 8 / 31 & 0.0416666666666667 $ pm $ 0.921569736993269 & NA & NA \\"
#[1] "100 & TQMC-n2 vs. ASTRAL3 & 10 / 8 / 30 & -0.0833333333333333 $ pm $ 1.48515108957481 & NA & NA \\"
#
#[1] "200 & TQMC-n2 vs. FASTRAL & 20 / 6 / 24 & 0.44 $ pm $ 1.09096474811639 & 6.216614e-03 & * \\"
#[1] "200 & TQMC-n2 vs. ASTRAL3 & 20 / 6 / 24 & 0.48 $ pm $ 1.12920416259254 & 4.216915e-03 & ** \\"
#
#[1] "500 & TQMC-n2 vs. FASTRAL & 34 / 9 / 7 & 2.4 $ pm $ 3.45820526768863 & 8.996877e-06 & ** \\"
#[1] "500 & TQMC-n2 vs. ASTRAL3 & 33 / 8 / 9 & 2.7 $ pm $ 3.89269331241106 & 1.269711e-05 & ** \\"
#
#[1] "1000 & TQMC-n2 vs. FASTRAL & 32 / 10 / 5 & 3.93617021276596 $ pm $ 6.95764907609244 & 1.227715e-04 & ** \\"
#[1] "1000 & TQMC-n2 vs. ASTRAL3 & 35 / 8 / 4 & 4.63829787234043 $ pm $ 7.83921450123032 & 7.717409e-06 & ** \\"
