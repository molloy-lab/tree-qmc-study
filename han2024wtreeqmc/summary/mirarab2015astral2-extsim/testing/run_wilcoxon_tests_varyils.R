#!/usr/bin/env Rscript

library(exactRankTests)
library(coin)

data <- read.csv("../csvs/data-varyils-for-testing.csv")

data$NTAX <- as.factor(data$NTAX)
data$ILSL <- as.factor(data$ILSL)
data$SPEC <- as.factor(data$SPEC)
data$SUPP <- as.factor(data$SUPP)
data$NGEN <- as.factor(data$NGEN)
data$REPL <- as.factor(data$REPL)

supps <- c("abayes")
ngens <- c(50, 200, 1000)

hghts <- c("500000", "2000000", "10000000")
ilsls <- c("low", "medium", "high")

rates <- c("1e-07", "1e-06")
specs <- c("deep", "shallow")

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
    message("on ASTRAL-II data for varying levels of ILS]{\\textbf{Testing for differences between TREE-QMC-wh\\_n2 vs")
    message(paste(mthd, "on the ASTRAL-II simulated data (with abayes support) for varying levels of ILS.} BET is the number of replicates for which"))
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
    message("\\begin{tabular}{r r r l r r r l c c c l}")
    message("\\toprule") 
    message(paste("\\multicolumn{12}{c}{\\textbf{TREE-QMC-wh\\_n2  vs ", mthd, "}} \\\\"))
    message("\\midrule")
    message("ILS Level & Speciation & \\# of genes & & BET & WOR & TIE & & p-val & sig & MC & note \\\\")
    message("\\midrule")

    for (supp in supps) {
        for (i in c(1, 2, 3)) {
                hght <- hghts[i]
                ilsl <- ilsls[i]
            for (j in c(1, 2)) {
                rate <- rates[j]
                spec <- specs[j]

                for (ngen in ngens) {

        #print(paste(as.character(i), as.character(j)))

        df <- data[data$SUPP == "abayes", ]  # FIXING
        df <- df[df$NGEN == ngen, ]
        df <- df[df$ILSL == ilsl, ]
        df <- df[df$SPEC == spec, ]

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

        writeme <- paste(ilsl, "&",
                         spec, "&",
                         as.character(ngen), "& &",
                         as.character(ntreeqmc), "&",
                         as.character(nother), "&",
                         as.character(ntie), "& &",
                         format(wsr, scientific = TRUE), "&",
                         "&", stars, "&", mc, note, "\\\\")
        message(writeme)

                }
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


\begin{table}[!h]
\caption[Testing for TREE-QMC-wh\_n2 vs CAML
on ASTRAL-II data for varying levels of ILS]{\textbf{Testing for differences between TREE-QMC-wh\_n2 vs
CAML on the ASTRAL-II simulated data (with abayes support) for varying levels of ILS.} BET is the number of replicates for which
TREE-QMC-wh\_n2 has lower species tree (RF) error and thus is better than CAML ,
WOR is the number of replicates for which wTREE-QMC-wh\_n2 has higher
RF error and thus is worse than CAML ,
and TIE is the number of replicates where the two methods tie.
Significance is evaluated using paired, two-sided Wilcoxon signed-rank tests on the RF error rates.
The symbols *, **, ***, ****, ***** indicate significance at $p <$ 0.5, 0.005, and so on.
MC indicates significance after Bonferroni correction, i.e., $p < $ 0.05 / 99 = 5e-04
for the  99 tests made on the S100 data.
}
\centering
\scriptsize
\begin{tabular}{r r r l r r r l c c c l}
\toprule
\multicolumn{12}{c}{\textbf{TREE-QMC-wh\_n2  vs  CAML }} \\
\midrule
ILS Level & Speciation & \# of genes & & BET & WOR & TIE & & p-val & sig & MC & note \\
\midrule
low & deep & 50 & & 12 & 31 & 7 & & 0.001014695 & & ** &  ( CAML better) \\
low & deep & 200 & & 13 & 28 & 9 & & 0.001628425 & & ** &  ( CAML better) \\
low & deep & 1000 & & 13 & 29 & 8 & & 0.01314316 & & * &  ( CAML better) \\
low & shallow & 50 & & 36 & 9 & 5 & & 5.632589e-06 & & **** & MC  \\
low & shallow & 200 & & 38 & 6 & 6 & & 8.63588e-06 & & **** & MC  \\
low & shallow & 1000 & & 29 & 10 & 11 & & 0.005499648 & & * &   \\
medium & deep & 50 & & 38 & 10 & 2 & & 4.729318e-06 & & ***** & MC  \\
medium & deep & 200 & & 37 & 9 & 4 & & 0.0004051359 & & *** & MC  \\
medium & deep & 1000 & & 24 & 19 & 7 & & 0.9784243 & &  &   \\
medium & shallow & 50 & & 44 & 3 & 3 & & 3.538503e-12 & & ***** & MC  \\
medium & shallow & 200 & & 40 & 4 & 6 & & 2.048751e-09 & & ***** & MC  \\
medium & shallow & 1000 & & 33 & 8 & 9 & & 0.0003452434 & & *** & MC  \\
high & deep & 50 & & 50 & 0 & 0 & & 1.776357e-15 & & ***** & MC  \\
high & deep & 200 & & 50 & 0 & 0 & & 1.776357e-15 & & ***** & MC  \\
high & deep & 1000 & & 45 & 4 & 0 & & 1.20962e-07 & & ***** & MC  \\
high & shallow & 50 & & 47 & 0 & 0 & & 1.421085e-14 & & ***** & MC  \\
high & shallow & 200 & & 47 & 0 & 0 & & 1.421085e-14 & & ***** & MC  \\
high & shallow & 1000 & & 46 & 0 & 1 & & 2.842171e-14 & & ***** & MC  \\
\bottomrule
\end{tabular}
\end{table}

\begin{table}[!h]
\caption[Testing for TREE-QMC-wh\_n2 vs ASTER-h
on ASTRAL-II data for varying levels of ILS]{\textbf{Testing for differences between TREE-QMC-wh\_n2 vs
ASTER-h on the ASTRAL-II simulated data (with abayes support) for varying levels of ILS.} BET is the number of replicates for which
TREE-QMC-wh\_n2 has lower species tree (RF) error and thus is better than ASTER-h ,
WOR is the number of replicates for which wTREE-QMC-wh\_n2 has higher
RF error and thus is worse than ASTER-h ,
and TIE is the number of replicates where the two methods tie.
Significance is evaluated using paired, two-sided Wilcoxon signed-rank tests on the RF error rates.
The symbols *, **, ***, ****, ***** indicate significance at $p <$ 0.5, 0.005, and so on.
MC indicates significance after Bonferroni correction, i.e., $p < $ 0.05 / 99 = 5e-04
for the  99 tests made on the S100 data.
}
\centering
\scriptsize
\begin{tabular}{r r r l r r r l c c c l}
\toprule
\multicolumn{12}{c}{\textbf{TREE-QMC-wh\_n2  vs  ASTER-h }} \\
\midrule
ILS Level & Speciation & \# of genes & & BET & WOR & TIE & & p-val & sig & MC & note \\
\midrule
low & deep & 50 & & 30 & 3 & 17 & & 2.060551e-07 & & ***** & MC  \\
low & deep & 200 & & 21 & 11 & 18 & & 0.02371259 & & * &   \\
low & deep & 1000 & & 26 & 2 & 22 & & 1.696497e-05 & & **** & MC  \\
low & shallow & 50 & & 32 & 10 & 8 & & 6.118376e-05 & & *** & MC  \\
low & shallow & 200 & & 21 & 13 & 16 & & 0.1606247 & &  &   \\
low & shallow & 1000 & & 17 & 5 & 28 & & 0.03724957 & & * &   \\
medium & deep & 50 & & 24 & 14 & 12 & & 0.05695275 & &  &   \\
medium & deep & 200 & & 29 & 8 & 13 & & 0.0001611809 & & *** & MC  \\
medium & deep & 1000 & & 22 & 11 & 17 & & 0.0186442 & & * &   \\
medium & shallow & 50 & & 22 & 17 & 11 & & 0.07322669 & &  &   \\
medium & shallow & 200 & & 25 & 10 & 15 & & 0.01056484 & & * &   \\
medium & shallow & 1000 & & 17 & 5 & 28 & & 0.01173639 & & * &   \\
high & deep & 50 & & 25 & 21 & 4 & & 0.2112505 & &  &   \\
high & deep & 200 & & 24 & 21 & 5 & & 0.9229223 & &  &   \\
high & deep & 1000 & & 20 & 17 & 12 & & 0.9970997 & &  &   \\
high & shallow & 50 & & 29 & 13 & 5 & & 0.01055602 & & * &   \\
high & shallow & 200 & & 29 & 12 & 6 & & 0.01537982 & & * &   \\
high & shallow & 1000 & & 26 & 14 & 7 & & 0.07923641 & &  &   \\
\bottomrule
\end{tabular}
\end{table}

\begin{table}[!h]
\caption[Testing for TREE-QMC-wh\_n2 vs wASTRID
on ASTRAL-II data for varying levels of ILS]{\textbf{Testing for differences between TREE-QMC-wh\_n2 vs
wASTRID on the ASTRAL-II simulated data (with abayes support) for varying levels of ILS.} BET is the number of replicates for which
TREE-QMC-wh\_n2 has lower species tree (RF) error and thus is better than wASTRID ,
WOR is the number of replicates for which wTREE-QMC-wh\_n2 has higher
RF error and thus is worse than wASTRID ,
and TIE is the number of replicates where the two methods tie.
Significance is evaluated using paired, two-sided Wilcoxon signed-rank tests on the RF error rates.
The symbols *, **, ***, ****, ***** indicate significance at $p <$ 0.5, 0.005, and so on.
MC indicates significance after Bonferroni correction, i.e., $p < $ 0.05 / 99 = 5e-04
for the  99 tests made on the S100 data.
}
\centering
\scriptsize
\begin{tabular}{r r r l r r r l c c c l}
\toprule
\multicolumn{12}{c}{\textbf{TREE-QMC-wh\_n2  vs  wASTRID }} \\
\midrule
ILS Level & Speciation & \# of genes & & BET & WOR & TIE & & p-val & sig & MC & note \\
\midrule
low & deep & 50 & & 31 & 7 & 12 & & 1.008374e-05 & & **** & MC  \\
low & deep & 200 & & 38 & 5 & 7 & & 2.443812e-09 & & ***** & MC  \\
low & deep & 1000 & & 43 & 2 & 5 & & 3.581135e-12 & & ***** & MC  \\
low & shallow & 50 & & 22 & 13 & 15 & & 0.2696331 & &  &   \\
low & shallow & 200 & & 12 & 20 & 18 & & 0.2902047 & &  &   \\
low & shallow & 1000 & & 14 & 7 & 29 & & 0.2642612 & &  &   \\
medium & deep & 50 & & 25 & 20 & 5 & & 0.05764566 & &  &   \\
medium & deep & 200 & & 33 & 10 & 7 & & 6.837434e-05 & & *** & MC  \\
medium & deep & 1000 & & 30 & 6 & 14 & & 4.754454e-05 & & **** & MC  \\
medium & shallow & 50 & & 30 & 13 & 7 & & 0.001086662 & & ** &   \\
medium & shallow & 200 & & 26 & 13 & 11 & & 0.02047486 & & * &   \\
medium & shallow & 1000 & & 26 & 12 & 12 & & 0.01566447 & & * &   \\
high & deep & 50 & & 35 & 9 & 6 & & 9.161746e-06 & & **** & MC  \\
high & deep & 200 & & 37 & 10 & 3 & & 1.476171e-05 & & **** & MC  \\
high & deep & 1000 & & 29 & 13 & 7 & & 0.1618871 & &  &   \\
high & shallow & 50 & & 39 & 5 & 3 & & 8.781171e-10 & & ***** & MC  \\
high & shallow & 200 & & 40 & 5 & 2 & & 5.238405e-08 & & ***** & MC  \\
high & shallow & 1000 & & 33 & 8 & 6 & & 1.492041e-06 & & ***** & MC  \\
\bottomrule
\end{tabular}
\end{table}
