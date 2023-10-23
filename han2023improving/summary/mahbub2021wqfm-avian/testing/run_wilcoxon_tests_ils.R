#!/usr/bin/env Rscript

# Remove ties
# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0200837

data <- read.csv("../csvs/data-for-testing.csv")
data$SCAL <- as.factor(data$SCAL)
data$NBPS <- as.factor(data$NBPS)
data$NGEN <- as.factor(data$NGEN)
data$REPL <- as.factor(data$REPL)

scals <- c("0.5X", "1X", "2X")
ngens <- c(1000)
nbpss <- c("500", "true")
mthds <- c("wQFM", "FASTRAL", "ASTRAL3")

total <- length(scals) * length(ngens) * length(nbpss) * length(mthds)

threshold2 <- 0.05
threshold1 <- 0.005
threshold0 <- threshold2 / total

writeme <- paste("* indicates p <",  as.character(threshold2))
print(writeme)
writeme <- paste("** indicates p <",  as.character(threshold1))
print(writeme)
writeme <- paste("Note: p <",  as.character(threshold0), "if correcting for", as.character(total), "comparisons")
print(writeme)

for (scal in scals) {
    xdf <- data[data$SCAL == scal, ]
    for (ngen in ngens) {
        qdf <- xdf[xdf$NGEN == ngen, ]
        for (nbps in nbpss) {
            ydf <- qdf[qdf$NBPS == nbps, ]
            for (mthd in mthds) {

                # Remove ties
                if (mthd == "wQFM") {
                    ydf$DIFF <- ydf$wQFMxSEFN - ydf$TQMCn2xSEFN  # positive means TQMCn2 is better
                } else if (mthd == "FASTRAL") {
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

                if (nrow(zdf) > 10) {
                    tmp1 <- zdf$TQMCn2xSEFN
                    if (mthd == "wQFM") {
                        tmp2 <- zdf$wQFMxSEFN
                    } else if (mthd == "FASTRAL") {
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
                        writeme <- paste(scal, "&",
                                         as.character(ngen), "&",
                                         as.character(nbps), "&",
                                         "TQMC-n2 vs.", mthd, "&", 
                                         as.character(ntreeqmc), "/",
                                         as.character(nother), "/",
                                         as.character(ntie), "&",
                                         as.character(avgdiffwtie), "$ pm $",
                                         as.character(vardiffwtie), "&",
                                         format(wsr$p.value, scientific = TRUE), "& ** \\")
                        print(writeme)
                    } else if (wsr$p.value < threshold2) {
                        writeme <- paste(scal, "&",
                                         as.character(ngen), "&",
                                         as.character(nbps), "&",
                                         "TQMC-n2 vs.", mthd, "&", 
                                         as.character(ntreeqmc), "/",
                                         as.character(nother), "/",
                                         as.character(ntie), "&",
                                         as.character(avgdiffwtie), "$ pm $",
                                         as.character(vardiffwtie), "&",
                                         format(wsr$p.value, scientific = TRUE), "& * \\")
                        print(writeme)
                    } else {
                        writeme <- paste(scal, "&",
                                         as.character(ngen), "&",
                                         as.character(nbps), "&",
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
                    writeme <- paste(scal, "&",
                                     as.character(ngen), "&",
                                     as.character(nbps), "&",
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
    }
}

#[1] "* indicates p < 0.05"
#[1] "** indicates p < 0.005"
#[1] "Note: p < 0.00277777777777778 if correcting for 18 comparisons"
#[1] "0.5X & 1000 & 500 & TQMC-n2 vs. wQFM & 11 / 2 / 7 & 0.75 $ pm $ 1.11803398874989 & 1.10058e-02 & * \\"
#[1] "0.5X & 1000 & 500 & TQMC-n2 vs. FASTRAL & 8 / 3 / 9 & 0.75 $ pm $ 1.48235232689554 & 4.605874e-02 & * \\"
#[1] "0.5X & 1000 & 500 & TQMC-n2 vs. ASTRAL3 & 18 / 2 / 0 & 2.1 $ pm $ 1.86095617754337 & 1.008441e-03 & ** \\"
#[1] "0.5X & 1000 & true & TQMC-n2 vs. wQFM & 5 / 7 / 8 & -0.05 $ pm $ 1.39453821823042 & 0.9359554 \\"
#[1] "0.5X & 1000 & true & TQMC-n2 vs. FASTRAL & 5 / 8 / 7 & -0.2 $ pm $ 1.23969435961577 & 0.5608834 \\"
#[1] "0.5X & 1000 & true & TQMC-n2 vs. ASTRAL3 & 9 / 7 / 4 & 0.15 $ pm $ 1.42441123571146 & 0.5570418 \\"
#[1] "1X & 1000 & 500 & TQMC-n2 vs. wQFM & 5 / 4 / 11 & 0.05 $ pm $ 0.68633274115326 & NA & NA \\"
#[1] "1X & 1000 & 500 & TQMC-n2 vs. FASTRAL & 9 / 3 / 8 & 0.7 $ pm $ 1.34164078649987 & 3.329306e-02 & * \\"
#[1] "1X & 1000 & 500 & TQMC-n2 vs. ASTRAL3 & 12 / 1 / 7 & 1.1 $ pm $ 1.2096106376586 & 2.925756e-03 & ** \\"
#[1] "1X & 1000 & true & TQMC-n2 vs. wQFM & 2 / 6 / 12 & -0.25 $ pm $ 0.716350399411379 & NA & NA \\"
#[1] "1X & 1000 & true & TQMC-n2 vs. FASTRAL & 2 / 7 / 11 & -0.3 $ pm $ 0.923380516876639 & NA & NA \\"
#[1] "1X & 1000 & true & TQMC-n2 vs. ASTRAL3 & 4 / 5 / 11 & -0.05 $ pm $ 0.887041208323017 & NA & NA \\"
#[1] "2X & 1000 & 500 & TQMC-n2 vs. wQFM & 4 / 5 / 11 & 0 $ pm $ 0.794719414239026 & NA & NA \\"
#[1] "2X & 1000 & 500 & TQMC-n2 vs. FASTRAL & 7 / 5 / 8 & 0.2 $ pm $ 0.951453182187509 & 0.3555315 \\"
#[1] "2X & 1000 & 500 & TQMC-n2 vs. ASTRAL3 & 7 / 5 / 8 & 0.15 $ pm $ 0.875093979915421 & 0.464435 \\"
#[1] "2X & 1000 & true & TQMC-n2 vs. wQFM & 0 / 0 / 20 & 0 $ pm $ 0 & NA & NA \\"
#[1] "2X & 1000 & true & TQMC-n2 vs. FASTRAL & 2 / 1 / 17 & 0.05 $ pm $ 0.394034462826206 & NA & NA \\"
#[1] "2X & 1000 & true & TQMC-n2 vs. ASTRAL3 & 2 / 1 / 17 & 0.05 $ pm $ 0.394034462826206 & NA & NA \\"
