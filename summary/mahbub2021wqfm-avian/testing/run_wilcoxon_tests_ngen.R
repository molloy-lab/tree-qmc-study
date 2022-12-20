#!/usr/bin/env Rscript

# Remove ties
# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0200837

data <- read.csv("../csvs/data-for-testing.csv")
data$SCAL <- as.factor(data$SCAL)
data$NBPS <- as.factor(data$NBPS)
data$NGEN <- as.factor(data$NGEN)
data$REPL <- as.factor(data$REPL)

scals <- c("1X")
ngens <- c(50, 100, 200, 500, 1000)
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


for (nbps in nbpss) {
    xdf <- data[data$NBPS == nbps, ]
    for (scal in scals) {
        qdf <- xdf[xdf$SCAL == scal, ]
        for (ngen in ngens) {
            ydf <- qdf[qdf$NGEN == ngen, ]
            
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
                        writeme <- paste(as.character(nbps), "&",
                                         scal, "&",
                                         as.character(ngen), "&",
                                         "TQMC-n2 vs.", mthd, "&", 
                                         as.character(ntreeqmc), "/",
                                         as.character(nother), "/",
                                         as.character(ntie), "&",
                                         as.character(avgdiffwtie), "$ pm $",
                                         as.character(vardiffwtie), "&",
                                         format(wsr$p.value, scientific = TRUE), "& ** \\")
                        print(writeme)
                    } else if (wsr$p.value < threshold2) {
                        writeme <- paste(as.character(nbps), "&",
                                         scal, "&",
                                         as.character(ngen), "&",
                                         "TQMC-n2 vs.", mthd, "&", 
                                         as.character(ntreeqmc), "/",
                                         as.character(nother), "/",
                                         as.character(ntie), "&",
                                         as.character(avgdiffwtie), "$ pm $",
                                         as.character(vardiffwtie), "&",
                                         format(wsr$p.value, scientific = TRUE), "& * \\")
                        print(writeme)
                    } else {
                        writeme <- paste(as.character(nbps), "&",
                                         scal, "&",
                                         as.character(ngen), "&",
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
                    writeme <- paste(as.character(nbps), "&",
                                     scal, "&",
                                     as.character(ngen), "&",
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
#[1] "Note: p < 0.00166666666666667 if correcting for 30 comparisons"

#[1] "500 & 1X & 50 & TQMC-n2 vs. wQFM & 3 / 8 / 9 & -0.35 $ pm $ 1.66306628661765 & 0.4716026 \\"
#[1] "500 & 1X & 50 & TQMC-n2 vs. FASTRAL & 11 / 5 / 4 & 0.8 $ pm $ 2.39736697671688 & 0.1578682 \\"
#[1] "500 & 1X & 50 & TQMC-n2 vs. ASTRAL3 & 8 / 6 / 6 & 0.35 $ pm $ 2.47673384244569 & 0.505146 \\"
#[1] "500 & 1X & 100 & TQMC-n2 vs. wQFM & 7 / 8 / 5 & 0.05 $ pm $ 1.66938375014948 & 0.9306406 \\"
#[1] "500 & 1X & 100 & TQMC-n2 vs. FASTRAL & 13 / 1 / 6 & 1.3 $ pm $ 1.49031964074119 & 2.023494e-03 & ** \\"
#[1] "500 & 1X & 100 & TQMC-n2 vs. ASTRAL3 & 13 / 2 / 5 & 1.1 $ pm $ 1.37266548230652 & 2.75909e-03 & ** \\"
#[1] "500 & 1X & 200 & TQMC-n2 vs. wQFM & 8 / 5 / 7 & 0.25 $ pm $ 1.40955386745706 & 0.3913659 \\"
#[1] "500 & 1X & 200 & TQMC-n2 vs. FASTRAL & 9 / 6 / 5 & 0.3 $ pm $ 1.52522647153585 & 0.3530914 \\"
#[1] "500 & 1X & 200 & TQMC-n2 vs. ASTRAL3 & 11 / 4 / 5 & 1 $ pm $ 1.74717817607346 & 2.66994e-02 & * \\"
#[1] "500 & 1X & 500 & TQMC-n2 vs. wQFM & 11 / 2 / 7 & 0.9 $ pm $ 1.4473205733718 & 1.93195e-02 & * \\"
#[1] "500 & 1X & 500 & TQMC-n2 vs. FASTRAL & 11 / 3 / 6 & 1.05 $ pm $ 1.8771478925557 & 2.632504e-02 & * \\"
#[1] "500 & 1X & 500 & TQMC-n2 vs. ASTRAL3 & 15 / 2 / 3 & 1.5 $ pm $ 1.50437957136384 & 1.120114e-03 & ** \\"
#[1] "500 & 1X & 1000 & TQMC-n2 vs. wQFM & 5 / 4 / 11 & 0.05 $ pm $ 0.68633274115326 & NA & NA \\"
#[1] "500 & 1X & 1000 & TQMC-n2 vs. FASTRAL & 9 / 3 / 8 & 0.7 $ pm $ 1.34164078649987 & 3.329306e-02 & * \\"
#[1] "500 & 1X & 1000 & TQMC-n2 vs. ASTRAL3 & 12 / 1 / 7 & 1.1 $ pm $ 1.2096106376586 & 2.925756e-03 & ** \\"

#[1] "true & 1X & 50 & TQMC-n2 vs. wQFM & 2 / 10 / 8 & -0.65 $ pm $ 1.18210338847862 & 3.714448e-02 & * \\"
#[1] "true & 1X & 50 & TQMC-n2 vs. FASTRAL & 11 / 5 / 4 & 0.7 $ pm $ 1.55935209217432 & 0.07447945 \\"
#[1] "true & 1X & 50 & TQMC-n2 vs. ASTRAL3 & 8 / 8 / 4 & 0 $ pm $ 1.4142135623731 & 1 \\"
#[1] "true & 1X & 100 & TQMC-n2 vs. wQFM & 7 / 6 / 7 & -0.2 $ pm $ 1.23969435961577 & 0.4672796 \\"
#[1] "true & 1X & 100 & TQMC-n2 vs. FASTRAL & 9 / 3 / 8 & 0.7 $ pm $ 1.55935209217432 & 0.07745305 \\"
#[1] "true & 1X & 100 & TQMC-n2 vs. ASTRAL3 & 9 / 4 / 7 & 0.6 $ pm $ 1.31389337066357 & 0.06917411 \\"
#[1] "true & 1X & 200 & TQMC-n2 vs. wQFM & 2 / 5 / 13 & -0.1 $ pm $ 0.967906041546987 & NA & NA \\"
#[1] "true & 1X & 200 & TQMC-n2 vs. FASTRAL & 8 / 5 / 7 & 0.2 $ pm $ 1.19648608323224 & 0.4926814 \\"
#[1] "true & 1X & 200 & TQMC-n2 vs. ASTRAL3 & 9 / 4 / 7 & 0.35 $ pm $ 1.22581873821025 & 0.2499018 \\"
#[1] "true & 1X & 500 & TQMC-n2 vs. wQFM & 5 / 0 / 15 & 0.3 $ pm $ 0.571240570577479 & NA & NA \\"
#[1] "true & 1X & 500 & TQMC-n2 vs. FASTRAL & 10 / 3 / 7 & 0.45 $ pm $ 0.887041208323017 & 4.2486e-02 & * \\"
#[1] "true & 1X & 500 & TQMC-n2 vs. ASTRAL3 & 10 / 3 / 7 & 0.45 $ pm $ 0.887041208323017 & 4.2486e-02 & * \\"
#[1] "true & 1X & 1000 & TQMC-n2 vs. wQFM & 2 / 6 / 12 & -0.25 $ pm $ 0.716350399411379 & NA & NA \\"
#[1] "true & 1X & 1000 & TQMC-n2 vs. FASTRAL & 2 / 7 / 11 & -0.3 $ pm $ 0.923380516876639 & NA & NA \\"
#[1] "true & 1X & 1000 & TQMC-n2 vs. ASTRAL3 & 4 / 5 / 11 & -0.05 $ pm $ 0.887041208323017 & NA & NA \\"
