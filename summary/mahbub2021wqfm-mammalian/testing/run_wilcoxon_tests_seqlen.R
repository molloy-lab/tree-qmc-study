#!/usr/bin/env Rscript

# Remove ties
# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0200837

data <- read.csv("../csvs/data-for-testing.csv")
data$SCAL <- as.factor(data$SCAL)
data$NBPS <- as.factor(data$NBPS)
data$NGEN <- as.factor(data$NGEN)
data$REPL <- as.factor(data$REPL)

scals <- c("noscale")
ngens <- c("200g")
nbpss <- c("250b", "500b", "1000b", "1500b", "true")
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
		        # 	print("Found a tie!")
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
#[1] "Note: p < 0.00333333333333333 if correcting for 15 comparisons"
#[1] "noscale & 200g & 250b & TQMC-n2 vs. wQFM & 0 / 0 / 20 & 0 $ pm $ 0 & NA & NA \\"
#[1] "noscale & 200g & 250b & TQMC-n2 vs. FASTRAL & 11 / 4 / 5 & 0.55 $ pm $ 1.09904264559757 & 4.303733e-02 & * \\"
#[1] "noscale & 200g & 250b & TQMC-n2 vs. ASTRAL3 & 11 / 4 / 5 & 0.55 $ pm $ 1.09904264559757 & 4.303733e-02 & * \\"
#[1] "noscale & 200g & 500b & TQMC-n2 vs. wQFM & 0 / 0 / 20 & 0 $ pm $ 0 & NA & NA \\"
#[1] "noscale & 200g & 500b & TQMC-n2 vs. FASTRAL & 5 / 2 / 13 & 0.15 $ pm $ 0.5871429486124 & NA & NA \\"
#[1] "noscale & 200g & 500b & TQMC-n2 vs. ASTRAL3 & 5 / 2 / 13 & 0.15 $ pm $ 0.5871429486124 & NA & NA \\"
#[1] "noscale & 200g & 1000b & TQMC-n2 vs. wQFM & 0 / 0 / 20 & 0 $ pm $ 0 & NA & NA \\"
#[1] "noscale & 200g & 1000b & TQMC-n2 vs. FASTRAL & 2 / 0 / 18 & 0.1 $ pm $ 0.307793505625546 & NA & NA \\"
#[1] "noscale & 200g & 1000b & TQMC-n2 vs. ASTRAL3 & 2 / 0 / 18 & 0.1 $ pm $ 0.307793505625546 & NA & NA \\"
#[1] "noscale & 200g & 1500b & TQMC-n2 vs. wQFM & 0 / 0 / 20 & 0 $ pm $ 0 & NA & NA \\"
#[1] "noscale & 200g & 1500b & TQMC-n2 vs. FASTRAL & 3 / 1 / 16 & 0.1 $ pm $ 0.447213595499958 & NA & NA \\"
#[1] "noscale & 200g & 1500b & TQMC-n2 vs. ASTRAL3 & 3 / 1 / 16 & 0.1 $ pm $ 0.447213595499958 & NA & NA \\"
#[1] "noscale & 200g & true & TQMC-n2 vs. wQFM & 0 / 0 / 20 & 0 $ pm $ 0 & NA & NA \\"
#[1] "noscale & 200g & true & TQMC-n2 vs. FASTRAL & 3 / 1 / 16 & 0.1 $ pm $ 0.447213595499958 & NA & NA \\"
#[1] "noscale & 200g & true & TQMC-n2 vs. ASTRAL3 & 3 / 1 / 16 & 0.1 $ pm $ 0.447213595499958 & NA & NA \\"
