#!/usr/bin/env Rscript

library(exactRankTests)
library(coin)

data <- read.csv("../csvs/data-varyntax-for-testing.csv")
data$NTAX <- as.factor(data$NTAX)
data$ILSL <- as.factor(data$ILSL)
data$SPEC <- as.factor(data$SPEC)
data$SUPP <- as.factor(data$SUPP)
data$NGEN <- as.factor(data$NGEN)
data$REPL <- as.factor(data$REPL)

supps <- c("abayes")
ngens <- c(50, 200, 1000)

ntaxs <- c(10, 50, 100, 200, 500, 1000)

mthds <- c("CAML", "ASTER-h", "wASTRID")

ntests <- length(supps) * length(ngens) * 11 * length(mthds)
threshold0 <- 0.05
threshold_bonferoni <- threshold0 / ntests

#writeme <- paste("alpha =", as.character(threshold0))
#print(writeme)
#print(paste("bonferroni alpha =", as.character(threshold_bonferoni), 
#            "=", as.character(threshold0), "/", as.character(ntests)))

for (mthd in mthds) {
    message("\\begin{table}[!h]")
    message(paste("\\caption[Testing for TREE-QMC-wh\\_n2 vs", mthd))
    message("on ASTRAL-II data for varying numbers of taxa]{\\textbf{Testing for differences between TREE-QMC-wh\\_n2 vs")
    message(paste(mthd, "on the ASTRAL-II simulated data (with abayes support) for varying numbers of taxa.} BET is the number of replicates for which"))
    message(paste("TREE-QMC-wh\\_n2 has lower species tree (RF) error and thus is better than", mthd, ",")) 
    message("WOR is the number of replicates for which wTREE-QMC-wh\\_n2 has higher")
    message(paste("RF error and thus is worse than", mthd, ",")) 
    message("and TIE is the number of replicates where the two methods tie.")
    message("Significance is evaluated using paired, two-sided Wilcoxon signed-rank tests on the RF error rates.")
    message("The symbols *, **, ***, ****, ***** indicate significance at $p <$ 0.5, 0.005, and so on.")
    message(paste("MC indicates significance after Bonferroni correction, i.e., $p < $",
                  as.character(threshold0), "/",
                  as.character(ntests), "=", 
                  format(threshold_bonferoni, scientific = TRUE, digits=1)))
    message(paste("for the ",
                  as.character(ntests),
                  "tests made on the S100 data."))
    message("}")
    message("\\centering")
    message("\\scriptsize")
    message("\\begin{tabular}{r r l r r r l c c c l}")
    message("\\toprule") 
    message(paste("\\multicolumn{11}{c}{\\textbf{TREE-QMC-wh\\_n2  vs ", mthd, "}} \\\\"))
    message("\\midrule")
    message("\\# of taxa & \\# of genes & & BET & WOR & TIE & & p-val & sig & MC & note \\\\")
    message("\\midrule")

    for (supp in supps) {
        for (ntax in ntaxs) {
            if ((mthd == "CAML") && (ntax == 1000)) {
                ngens <- c(50, 200)
            } else {
                ngens <- c(50, 200, 1000)
            }
            for (ngen in ngens) {

        #writeme <- paste(mthd, supp, as.character(ntax), as.character(ngen))
        #print(writeme)

        fndf <- data[data$SUPP == supp, ]
        fndf <- fndf[fndf$NGEN == ngen, ]
        fndf <- fndf[fndf$NTAX == ntax, ]

        # Find ties
        mthd1 <-fndf$TQMCwhn2xSEFN
        if (mthd == "CAML") {
            mthd2 <- fndf$CAMLxSEFN
        } else if (mthd == "ASTER-h") {
            mthd2 <- fndf$ASTERHxSEFN
        } else if (mthd == "wASTRID") {
            mthd2 <- fndf$WASTRIDxSEFN
        } else {
            print("ERROR")
            exit()
        }

		#if (0 %in% fndf$DIFF) {
		#	print("Found a tie!")
		#}

        # Count ties
        diff <- mthd2 - mthd1
        ntreeqmc <- sum(diff > 0)
		nother <- sum(diff < 0)
        ntie <- sum(diff == 0)

        # Compute average difference
        #avgdiffwtie <- mean(ydf$DIFF)
        #vardiffwtie <- sd(ydf$DIFF)
        #meddiffwtie <- median(ydf$DIFF)

        #avgdiff <- mean(zdf$DIFF)
        #meddiff <- median(zdf$DIFF)

        # Run test
        if (ntreeqmc + nother == 0) {
            writeme <- paste(as.character(ntax), "&",
                             as.character(ngen), "& &",
                             as.character(ntreeqmc), "&",
                             as.character(nother), "&",
                             as.character(ntie), "& &",
                             "NA & NA & NA \\\\")
        } else {
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
            writeme <- paste(as.character(ntax), "&",
                             as.character(ngen), "& &",
                             as.character(ntreeqmc), "&",
                             as.character(nother), "&",
                             as.character(ntie), "& &",
                             format(wsr, scientific = TRUE), "&", 
                             stars, "&", mc, note, "\\\\")

        }
        message(writeme)

            }
        }
    }
    message("\\bottomrule")
    message("\\end{tabular}")
    message("\\end{table}")
    message("")
}

print("done")
exit()

[1] "alpha = 0.05"
[1] "bonferroni alpha = 0.000505050505050505 = 0.05 / 99"
[1] "TQMC-wh-n2 vs. CAML & abayes & 10 & 50 & 8 / 6 / 35 & 0.6334229 & 0.6334229 & &  &  "
[1] "TQMC-wh-n2 vs. CAML & abayes & 10 & 200 & 4 / 2 / 43 & 0.53125 & 0.53125 & &  &  "
[1] "TQMC-wh-n2 vs. CAML & abayes & 10 & 1000 & 2 / 1 / 45 & 0.75 & 0.75 & &  &  "
[1] "TQMC-wh-n2 vs. CAML & abayes & 50 & 50 & 28 / 9 / 10 & 0.0001215681 & 0.0001925013 & & *** & MC "
[1] "TQMC-wh-n2 vs. CAML & abayes & 50 & 200 & 29 / 12 / 6 & 0.0006905201 & 0.0009076932 & & ** &  "
[1] "TQMC-wh-n2 vs. CAML & abayes & 50 & 1000 & 21 / 3 / 23 & 0.0007148981 & 0.0002533197 & & ** &  "
[1] "TQMC-wh-n2 vs. CAML & abayes & 100 & 50 & 38 / 5 / 5 & 1.594049e-08 & 1.284025e-08 & & ***** & MC "
[1] "TQMC-wh-n2 vs. CAML & abayes & 100 & 200 & 32 / 9 / 7 & 0.001508224 & 0.0008027762 & & ** &  "
[1] "TQMC-wh-n2 vs. CAML & abayes & 100 & 1000 & 23 / 13 / 12 & 0.03728368 & 0.04959638 & & * &  "
[1] "TQMC-wh-n2 vs. CAML & abayes & 200 & 50 & 44 / 3 / 3 & 3.012701e-12 & 3.538503e-12 & & ***** & MC "
[1] "TQMC-wh-n2 vs. CAML & abayes & 200 & 200 & 40 / 4 / 6 & 2.048751e-09 & 1.34969e-09 & & ***** & MC "
[1] "TQMC-wh-n2 vs. CAML & abayes & 200 & 1000 & 33 / 8 / 9 & 0.0003452434 & 0.0001578385 & & *** & MC "
[1] "TQMC-wh-n2 vs. CAML & abayes & 500 & 50 & 46 / 3 / 1 & 6.039613e-14 & 8.881784e-14 & & ***** & MC "
[1] "TQMC-wh-n2 vs. CAML & abayes & 500 & 200 & 42 / 7 / 1 & 3.823587e-09 & 3.769735e-09 & & ***** & MC "
[1] "TQMC-wh-n2 vs. CAML & abayes & 500 & 1000 & 36 / 11 / 3 & 0.004434673 & 0.003003294 & & ** &  "
[1] "TQMC-wh-n2 vs. CAML & abayes & 1000 & 50 & 45 / 5 / 0 & 8.173611e-08 & 8.173611e-08 & & ***** & MC "
[1] "TQMC-wh-n2 vs. CAML & abayes & 1000 & 200 & 42 / 8 / 0 & 1.173318e-06 & 1.173318e-06 & & ***** & MC "
[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 10 & 50 & 1 / 1 / 47 & 1 & 1 & &  &  "
[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 10 & 200 & 0 / 1 / 48 & 1 & 1 & &  &  "
[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 10 & 1000 & 0 / 0 / 48 & NA & NA"
[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 50 & 50 & 11 / 8 / 28 & 0.324913 & 0.4520187 & &  &  "
[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 50 & 200 & 8 / 3 / 36 & 0.2265625 & 0.2265625 & &  &  "
[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 50 & 1000 & 5 / 3 / 39 & 0.7265625 & 0.7265625 & &  &  "
[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 100 & 50 & 17 / 8 / 23 & 0.01875782 & 0.0395084 & & * &  "
[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 100 & 200 & 15 / 7 / 26 & 0.1132421 & 0.1030617 & &  &  "
[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 100 & 1000 & 8 / 10 / 30 & 0.6668396 & 0.67379 & &  &  "
[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 200 & 50 & 22 / 17 / 11 & 0.02650681 & 0.07322669 & &  &  "
[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 200 & 200 & 25 / 10 / 15 & 0.01056484 & 0.00998523 & & * &  "
[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 200 & 1000 & 17 / 5 / 28 & 0.01054907 & 0.01173639 & & * &  "
[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 500 & 50 & 38 / 10 / 2 & 1.453705e-07 & 1.885201e-07 & & ***** & MC "
[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 500 & 200 & 40 / 6 / 4 & 4.758704e-09 & 5.277997e-09 & & ***** & MC "
[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 500 & 1000 & 36 / 8 / 6 & 1.304655e-05 & 8.947205e-06 & & **** & MC "
[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 1000 & 50 & 44 / 2 / 4 & 3.126388e-13 & 4.547474e-13 & & ***** & MC "
[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 1000 & 200 & 38 / 4 / 8 & 1.718945e-10 & 4.101821e-10 & & ***** & MC "
[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 1000 & 1000 & 44 / 2 / 4 & 5.684342e-12 & 3.524292e-12 & & ***** & MC "
[1] "TQMC-wh-n2 vs. wASTRID & abayes & 10 & 50 & 1 / 2 / 46 & 1 & 1 & &  &  "
[1] "TQMC-wh-n2 vs. wASTRID & abayes & 10 & 200 & 0 / 1 / 48 & 1 & 1 & &  &  "
[1] "TQMC-wh-n2 vs. wASTRID & abayes & 10 & 1000 & 0 / 0 / 48 & NA & NA"
[1] "TQMC-wh-n2 vs. wASTRID & abayes & 50 & 50 & 12 / 13 / 22 & 0.8214659 & 0.8604026 & &  &  "
[1] "TQMC-wh-n2 vs. wASTRID & abayes & 50 & 200 & 18 / 3 / 26 & 0.00435257 & 0.001437187 & & ** &  "
[1] "TQMC-wh-n2 vs. wASTRID & abayes & 50 & 1000 & 6 / 3 / 38 & 0.5078125 & 0.5078125 & &  &  "
[1] "TQMC-wh-n2 vs. wASTRID & abayes & 100 & 50 & 21 / 14 / 13 & 0.1207092 & 0.1545381 & &  &  "
[1] "TQMC-wh-n2 vs. wASTRID & abayes & 100 & 200 & 19 / 13 / 16 & 0.1220586 & 0.1799249 & &  &  "
[1] "TQMC-wh-n2 vs. wASTRID & abayes & 100 & 1000 & 14 / 7 / 27 & 0.02131462 & 0.0765667 & &  &  "
[1] "TQMC-wh-n2 vs. wASTRID & abayes & 200 & 50 & 30 / 13 / 7 & 0.0007217834 & 0.001086662 & & ** &  "
[1] "TQMC-wh-n2 vs. wASTRID & abayes & 200 & 200 & 26 / 13 / 11 & 0.01845867 & 0.02047486 & & * &  "
[1] "TQMC-wh-n2 vs. wASTRID & abayes & 200 & 1000 & 26 / 12 / 12 & 0.01553998 & 0.01566447 & & * &  "
[1] "TQMC-wh-n2 vs. wASTRID & abayes & 500 & 50 & 35 / 13 / 2 & 0.0001342302 & 0.0001402717 & & *** & MC "
[1] "TQMC-wh-n2 vs. wASTRID & abayes & 500 & 200 & 41 / 7 / 2 & 4.245093e-08 & 3.9225e-08 & & ***** & MC "
[1] "TQMC-wh-n2 vs. wASTRID & abayes & 500 & 1000 & 29 / 9 / 12 & 8.476483e-06 & 3.824975e-05 & & **** & MC "
[1] "TQMC-wh-n2 vs. wASTRID & abayes & 1000 & 50 & 43 / 6 / 1 & 1.760259e-09 & 1.654236e-09 & & ***** & MC "
[1] "TQMC-wh-n2 vs. wASTRID & abayes & 1000 & 200 & 46 / 3 / 1 & 7.94742e-12 & 6.000533e-12 & & ***** & MC "
[1] "TQMC-wh-n2 vs. wASTRID & abayes & 1000 & 1000 & 41 / 6 / 3 & 9.718804e-11 & 1.525962e-10 & & ***** & MC "
[1] "done"
