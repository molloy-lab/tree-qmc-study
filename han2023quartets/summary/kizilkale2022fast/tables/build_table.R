#!/usr/bin/env Rscript

# Remove ties
# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0200837

data <- read.csv("../csvs/data-for-table.csv")
data$NCELL <- as.factor(data$NCELL)
data$NMUT <- as.factor(data$NMUT)
data$REPL <- as.factor(data$REPL)

conds <- c(1, 2, 3)
ncells <- c("n1000", "n300", "n300")
nmuts <- c("m300", "m300", "m1000")
mthds <- c("FASTRAL", "SCISTREE", "FASTME")

total <- length(conds)

threshold2 <- 0.05
threshold1 <- 0.005
threshold0 <- threshold2 / total

writeme <- paste("* indicates p <",  as.character(threshold2))
print(writeme)
writeme <- paste("** indicates p <",  as.character(threshold1))
print(writeme)
writeme <- paste("Note: p <",  as.character(threshold0), "if correcting for", as.character(total), "comparisons")
print(writeme)


for (cond in conds) {
    ncell <- ncells[cond]
    nmut <- nmuts[cond]

    xdf <- data[data$NCELL == ncell, ]
    ydf <- xdf[xdf$NMUT == nmut, ]
    for (mthd in mthds) {
        # Remove ties
        if (mthd == "FASTRAL") {
            ydf$DIFF <- ydf$TQMCn2xNQS - ydf$FASTRALxNQS  # positive means TQMCn2 is better
        } else if (mthd == "SCISTREE") {
            ydf$DIFF <- ydf$TQMCn2xNQS - ydf$SCISTREExNQS # positive means TQMCn2 is better
        } else {
            ydf$DIFF <- ydf$TQMCn2xNQS - ydf$FASTMExNQS   # positive means TQMCn2 is better
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

        if (nrow(zdf) == 10) {
            tmp1 <- zdf$TQMCn2xNQS
            if (mthd == "FASTRAL") {
                tmp2 <- zdf$FASTRALxNQS
            } else if (mthd == "SCISTREE") {
                tmp2 <- zdf$SCISTREExNQS
            } else {
                tmp2 <- zdf$FASTMExNQS
            }

            wsr <- wilcox.test(tmp1, 
                               tmp2,
                               paired=TRUE,
                               alternative="two.sided",
                               mu=0, exact=NULL)
            #print(wsr)

            if (wsr$p.value < threshold1) {
                writeme <- paste(as.character(ncell), "&",
                                 as.character(nmut), "&",
                                 "TQMC-n2 vs.", mthd, "&", 
                                 as.character(ntreeqmc), "/",
                                 as.character(nother), "/",
                                 as.character(ntie), "&",
                                 as.character(avgdiffwtie), "$ pm $",
                                 as.character(vardiffwtie), "&",
                                 format(wsr$p.value, scientific = TRUE), "& ** \\")
                print(writeme)
             } else if (wsr$p.value < threshold2) {
                writeme <- paste(as.character(ncell), "&",
                                 as.character(nmut), "&",
                                 "TQMC-n2 vs.", mthd, "&", 
                                 as.character(ntreeqmc), "/",
                                 as.character(nother), "/",
                                 as.character(ntie), "&",
                                 as.character(avgdiffwtie), "$ pm $",
                                 as.character(vardiffwtie), "&",
                                 format(wsr$p.value, scientific = TRUE), "& * \\")
                print(writeme)
            } else {
                writeme <- paste(as.character(ncell), "&",
                                 as.character(nmut), "&",
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
            writeme <- paste(as.character(ncell), "&",
                             as.character(nmut), "&",
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
#[1] "Note: p < 0.0166666666666667 if correcting for 3 comparisons"
#[1] "n1000 & m300 & TQMC-n2 vs. FASTRAL & 10 / 0 / 0 & 0.00360149684972026 $ pm $ 0.00336157492387845 & 1.953125e-03 & ** \\"
#[1] "n1000 & m300 & TQMC-n2 vs. SCISTREE & 10 / 0 / 0 & 0.00622365945527301 $ pm $ 0.00833996044205936 & 1.953125e-03 & ** \\"
#[1] "n1000 & m300 & TQMC-n2 vs. FASTME & 10 / 0 / 0 & 0.00361197938370648 $ pm $ 0.00393847668144443 & 1.953125e-03 & ** \\"
#[1] "n300 & m300 & TQMC-n2 vs. FASTRAL & 5 / 5 / 0 & 0.000791133990277804 $ pm $ 0.00292555244216231 & 0.4316406 \\"
#[1] "n300 & m300 & TQMC-n2 vs. SCISTREE & 7 / 3 / 0 & 0.00121691994550456 $ pm $ 0.00450281649071144 & 0.3222656 \\"
#[1] "n300 & m300 & TQMC-n2 vs. FASTME & 8 / 2 / 0 & 0.00228831273133417 $ pm $ 0.00388638650138969 & 0.06445312 \\"
#[1] "n300 & m1000 & TQMC-n2 vs. FASTRAL & 0 / 10 / 0 & -0.000460470477324093 $ pm $ 0.000182239199116249 & 1.953125e-03 & ** \\"
#[1] "n300 & m1000 & TQMC-n2 vs. SCISTREE & 6 / 4 / 0 & 0.00031143579216566 $ pm $ 0.000495846962478827 & 0.1601562 \\"
#[1] "n300 & m1000 & TQMC-n2 vs. FASTME & 9 / 1 / 0 & 0.000758769441831075 $ pm $ 0.00061756447374818 & 3.90625e-03 & ** \\"
