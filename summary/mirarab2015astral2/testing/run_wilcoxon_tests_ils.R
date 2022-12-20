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

data <- subset(data, NTAX == 200)
data <- subset(data, GTRE == "estimatedgenetre")
data <- subset(data, NGEN == 1000)

strhts <- c("500000", "2000000", "10000000")
ilss <- c("0.25X", "1X", "5X")
srates <- c("1e-07", "1e-06")
speciations <- c("deep", "shallow")
mthds <- c("FASTRAL", "ASTRAL3")

total <- length(strhts) * length(srates) * length(mthds)

threshold2 <- 0.05
threshold1 <- 0.005
threshold0 <- threshold2 / total

writeme <- paste("* indicates p <",  as.character(threshold2))
print(writeme)
writeme <- paste("** indicates p <",  as.character(threshold1))
print(writeme)
writeme <- paste("Note: p <",  as.character(threshold0), "if correcting for", as.character(total), "comparisons")
print(writeme)


for (i in c(1, 2, 3)) {
    strht <- strhts[i]
    ils <- ilss[i]
    xdf <- data[data$STRHT == strht, ]
    
    for (j in c(1, 2)) {
        #print(paste(as.character(i), as.character(j)))
        srate <- srates[j]
        speciation <- speciations[j]
        ydf <- xdf[xdf$SRATE == srate, ]

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
                    writeme <- paste(ils, "&",
                                     speciation, "&",
                                     "TQMC-n2 vs.", mthd, "&", 
                                     as.character(ntreeqmc), "/",
                                     as.character(nother), "/",
                                     as.character(ntie), "&",
                                     as.character(avgdiffwtie), "$ pm $",
                                     as.character(vardiffwtie), "&",
                                     format(wsr$p.value, scientific = TRUE), "& ** \\")
                    print(writeme)
                } else if (wsr$p.value < threshold2) {
                    writeme <- paste(ils, "&",
                                     speciation, "&",
                                     "TQMC-n2 vs.", mthd, "&", 
                                     as.character(ntreeqmc), "/",
                                     as.character(nother), "/",
                                     as.character(ntie), "&",
                                     as.character(avgdiffwtie), "$ pm $",
                                     as.character(vardiffwtie), "&",
                                     format(wsr$p.value, scientific = TRUE), "& * \\")
                    print(writeme)
                } else {
                    writeme <- paste(ils, "&",
                                     speciation, "&",
                                     "TQMC-n2 vs.", mthd, "&", 
                                     as.character(ntreeqmc), "/",
                                     as.character(nother), "/",
                                     as.character(ntie), "&",
                                     as.character(avgdiffwtie), "$ pm $",
                                     as.character(vardiffwtie), "&",
                                     format(wsr$p.value, scientific = FALSE), "& \\")
                    print(writeme)
                }
            } else {
                writeme <- paste(ils, "&",
                                 speciation, "&",
                                 "TQMC-n2 vs.", mthd, "&", 
                                 as.character(ntreeqmc), "/",
                                 as.character(nother), "/",
                                 as.character(ntie), "&",
                                 as.character(avgdiffwtie), "$ pm $",
                                 as.character(vardiffwtie),
                                 "& NA & NA \\")
                print(writeme)
            }
        }
    }
}

#[1] "* indicates p < 0.05"
#[1] "** indicates p < 0.005"
#[1] "Note: p < 0.00416666666666667 if correcting for 12 comparisons"
#[1] "0.25X & deep & TQMC-n2 vs. FASTRAL & 22 / 21 / 7 & -0.24 $ pm $ 3.46150845966416 & 0.701991 & \\"
#[1] "0.25X & deep & TQMC-n2 vs. ASTRAL3 & 29 / 18 / 3 & 0.5 $ pm $ 3.52397039529829 & 0.09590108 & \\"
#[1] "0.25X & shallow & TQMC-n2 vs. FASTRAL & 23 / 16 / 8 & 0.468085106382979 $ pm $ 2.24432679777437 & 0.1559253 & \\"
#[1] "0.25X & shallow & TQMC-n2 vs. ASTRAL3 & 23 / 15 / 9 & 0.765957446808511 $ pm $ 2.3518056518847 & 4.41384e-02 & * \\"
#[1] "1X & deep & TQMC-n2 vs. FASTRAL & 20 / 9 / 21 & 0.7 $ pm $ 2.14998813475606 & 4.20894e-02 & * \\"
#[1] "1X & deep & TQMC-n2 vs. ASTRAL3 & 21 / 7 / 22 & 1 $ pm $ 2.24971653543195 & 3.666061e-03 & ** \\"
#[1] "1X & shallow & TQMC-n2 vs. FASTRAL & 20 / 6 / 24 & 0.44 $ pm $ 1.09096474811639 & 6.216614e-03 & * \\"
#[1] "1X & shallow & TQMC-n2 vs. ASTRAL3 & 20 / 6 / 24 & 0.48 $ pm $ 1.12920416259254 & 4.216915e-03 & ** \\"
#[1] "5X & deep & TQMC-n2 vs. FASTRAL & 27 / 6 / 17 & 1.1 $ pm $ 1.87627507665196 & 3.165028e-04 & ** \\"
#[1] "5X & deep & TQMC-n2 vs. ASTRAL3 & 25 / 8 / 17 & 0.9 $ pm $ 2.26102739876349 & 3.818722e-03 & ** \\"
#[1] "5X & shallow & TQMC-n2 vs. FASTRAL & 13 / 8 / 29 & 0.44 $ pm $ 1.38740691910617 & NA & NA \\"
#[1] "5X & shallow & TQMC-n2 vs. ASTRAL3 & 13 / 9 / 28 & 0.42 $ pm $ 1.41551173410458 & NA & NA \\"
