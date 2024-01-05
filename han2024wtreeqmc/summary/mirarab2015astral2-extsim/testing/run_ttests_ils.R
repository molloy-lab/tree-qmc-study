#!/usr/bin/env Rscript

# Remove ties
# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0200837



data <- read.csv("../csvs/data-for-testing.csv")

data$NTAX <- as.factor(data$NTAX)
data$HGHT <- as.factor(data$HGHT)
data$RATE <- as.factor(data$RATE)
data$SUPP <- as.factor(data$SUPP)
data$NGEN <- as.factor(data$NGEN)
data$REPL <- as.factor(data$REPL)

data <- subset(data, NTAX == 200)

supps <- c("abayes")
ngens <- c(50, 200, 1000)

hghts <- c("500000", "2000000", "10000000")
ilss <- c("0.25X", "1X", "5X")

rates <- c("1e-07", "1e-06")
specs <- c("deep", "shallow")

mthds <- c("CAML", "ASTER-h", "wASTRID")

ntest <- length(supps) * length(ngens) * 
         length(hghts) * length(rates) * 
         length(mthds)
threshold_bonferoni <- 0.05 / ntest

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

        fndf <- data[data$SUPP == supp, ]
        fndf <- fndf[fndf$NGEN == ngen, ]
        fndf <- fndf[fndf$HGHT == hght, ]
        fndf <- fndf[fndf$RATE == rate, ]

        # Find ties
        if (mthd == "CAML") {
            fndf$DIFF <- fndf$CAMLxSEFN - fndf$TQMCwhn2xSEFN  # positive means TQMC-wh-n2 is better
        } else if (mthd == "ASTER-h") {
            fndf$DIFF <- fndf$ASTERHxSEFN - fndf$TQMCwhn2xSEFN  # positive means TQMC-wh-n2 is better
        } else if (mthd == "wASTRID") {
            fndf$DIFF <- fndf$WASTRIDxSEFN - fndf$TQMCwhn2xSEFN  # positive means TQMC-wh-n2 is better
        } else {
            print("ERROR")
            exit()
        }

        # Count ties
        ntie <- sum(fndf$DIFF == 0)
        nother <- sum(fndf$DIFF < 0)
        ntreeqmc <- sum(fndf$DIFF > 0)

        # Run test
        tmp1 <- fndf$TQMCwhn2xSEFN
        if (mthd == "CAML") {
            tmp2 <- fndf$CAMLxSEFN
        } else if (mthd == "ASTER-h") {
            tmp2 <- fndf$ASTERHxSEFN
        } else if (mthd == "wASTRID") {
            tmp2 <- fndf$WASTRIDxSEFN
        } else {
            print("ERROR")
            exit()
        }

        tt <- t.test(tmp1, y = tmp2,
                     alternative="two.sided",
                     mu = 0, paired=TRUE, var.equal=FALSE,
                     conf.level=0.95)

        # Report
        if (ntreeqmc + nother == 0) {
            # ONLY TIESSSSS
            tt$p.value <- "NA"
            stars = ""
            mc <- ""
        } else {
            if (tt$p.value < 0.000005) {
                stars = "*****"
            } else if (tt$p.value < 0.00005) {
                stars = "****"
            } else if (tt$p.value < 0.0005) {
                stars = "***"
            } else if (tt$p.value < 0.005) {
                stars = "**"
            } else if (tt$p.value < 0.05) {
                stars = "*"
            } else {
                stars = ""
            }
            if (tt$p.value < threshold_bonferoni) {
                mc <- "MC"
            } else {
                mc <- ""
            }
        }

        writeme <- paste("TQMC-wh-n2 vs.", mthd, "&",
                         supp, "&",
                         ils, "&",
                         spec, "&",
                         as.character(ngen), "&",
                         as.character(ntreeqmc), "/",
                         as.character(nother), "/",
                         as.character(ntie), "&",
                         format(tt$p.value, scientific = TRUE), "&", stars, "&", mc)
        print(writeme)

                }
            }
        }
    }
}

#[1] "TQMC-wh-n2 vs. CAML & abayes & 0.25X & deep & 50 & 50 / 0 / 0 & 3.248776e-22 & ***** & MC"
#[1] "TQMC-wh-n2 vs. CAML & abayes & 0.25X & deep & 200 & 50 / 0 / 0 & 6.534075e-18 & ***** & MC"
#[1] "TQMC-wh-n2 vs. CAML & abayes & 0.25X & deep & 1000 & 46 / 4 / 0 & 3.137395e-08 & ***** & MC"
#[1] "TQMC-wh-n2 vs. CAML & abayes & 0.25X & shallow & 50 & 47 / 0 / 0 & 1.305157e-23 & ***** & MC"
#[1] "TQMC-wh-n2 vs. CAML & abayes & 0.25X & shallow & 200 & 47 / 0 / 0 & 9.466453e-16 & ***** & MC"
#[1] "TQMC-wh-n2 vs. CAML & abayes & 0.25X & shallow & 1000 & 46 / 0 / 1 & 1.099124e-16 & ***** & MC"
#[1] "TQMC-wh-n2 vs. CAML & abayes & 1X & deep & 50 & 38 / 10 / 2 & 8.120213e-05 & *** & MC"
#[1] "TQMC-wh-n2 vs. CAML & abayes & 1X & deep & 200 & 37 / 9 / 4 & 6.159305e-03 & * & "
#[1] "TQMC-wh-n2 vs. CAML & abayes & 1X & deep & 1000 & 24 / 19 / 7 & 4.045701e-01 &  & "
#[1] "TQMC-wh-n2 vs. CAML & abayes & 1X & shallow & 50 & 44 / 3 / 3 & 5.084971e-12 & ***** & MC"
#[1] "TQMC-wh-n2 vs. CAML & abayes & 1X & shallow & 200 & 40 / 4 / 6 & 3.090232e-09 & ***** & MC"
#[1] "TQMC-wh-n2 vs. CAML & abayes & 1X & shallow & 1000 & 33 / 8 / 9 & 1.054732e-03 & ** & "
#[1] "TQMC-wh-n2 vs. CAML & abayes & 5X & deep & 50 & 12 / 31 / 7 & 5.88925e-04 & ** & MC"
#[1] "TQMC-wh-n2 vs. CAML & abayes & 5X & deep & 200 & 13 / 28 / 9 & 1.404269e-03 & ** & "
#[1] "TQMC-wh-n2 vs. CAML & abayes & 5X & deep & 1000 & 13 / 29 / 8 & 2.623325e-02 & * & "
#[1] "TQMC-wh-n2 vs. CAML & abayes & 5X & shallow & 50 & 36 / 9 / 5 & 3.88531e-06 & ***** & MC"
#[1] "TQMC-wh-n2 vs. CAML & abayes & 5X & shallow & 200 & 38 / 6 / 6 & 1.673949e-05 & **** & MC"
#[1] "TQMC-wh-n2 vs. CAML & abayes & 5X & shallow & 1000 & 29 / 10 / 11 & 2.421272e-02 & * & "

#[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 0.25X & deep & 50 & 25 / 21 / 4 & 1.992515e-01 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 0.25X & deep & 200 & 24 / 21 / 5 & 7.927034e-01 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 0.25X & deep & 1000 & 21 / 17 / 12 & 5.374779e-01 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 0.25X & shallow & 50 & 29 / 13 / 5 & 2.612718e-02 & * & "
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 0.25X & shallow & 200 & 29 / 12 / 6 & 5.728345e-03 & * & "
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 0.25X & shallow & 1000 & 26 / 14 / 7 & 8.319004e-02 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 1X & deep & 50 & 24 / 14 / 12 & 3.610141e-02 & * & "
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 1X & deep & 200 & 29 / 8 / 13 & 1.814297e-04 & *** & MC"
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 1X & deep & 1000 & 22 / 11 / 17 & 7.073938e-03 & * & "
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 1X & shallow & 50 & 22 / 17 / 11 & 1.697341e-02 & * & "
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 1X & shallow & 200 & 25 / 10 / 15 & 1.244558e-02 & * & "
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 1X & shallow & 1000 & 17 / 5 / 28 & 6.8334e-03 & * & "
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 5X & deep & 50 & 30 / 3 / 17 & 8.779537e-07 & ***** & MC"
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 5X & deep & 200 & 21 / 11 / 18 & 1.424496e-02 & * & "
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 5X & deep & 1000 & 26 / 2 / 22 & 3.3288e-04 & *** & MC"
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 5X & shallow & 50 & 32 / 10 / 8 & 5.982737e-05 & *** & MC"
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 5X & shallow & 200 & 21 / 13 / 16 & 1.39581e-01 &  & "

#[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 5X & shallow & 1000 & 17 / 5 / 28 & 2.53029e-02 & * & "
#[1] "TQMC-wh-n2 vs. wASTRID & abayes & 0.25X & deep & 50 & 35 / 9 / 6 & 1.303658e-05 & **** & MC"
#[1] "TQMC-wh-n2 vs. wASTRID & abayes & 0.25X & deep & 200 & 37 / 10 / 3 & 8.925647e-06 & **** & MC"
#[1] "TQMC-wh-n2 vs. wASTRID & abayes & 0.25X & deep & 1000 & 30 / 13 / 7 & 7.237602e-01 &  & "
#[1] "TQMC-wh-n2 vs. wASTRID & abayes & 0.25X & shallow & 50 & 39 / 5 / 3 & 4.927677e-10 & ***** & MC"
#[1] "TQMC-wh-n2 vs. wASTRID & abayes & 0.25X & shallow & 200 & 40 / 5 / 2 & 1.748761e-08 & ***** & MC"
#[1] "TQMC-wh-n2 vs. wASTRID & abayes & 0.25X & shallow & 1000 & 33 / 8 / 6 & 5.702898e-07 & ***** & MC"
#[1] "TQMC-wh-n2 vs. wASTRID & abayes & 1X & deep & 50 & 25 / 20 / 5 & 1.608211e-02 & * & "
#[1] "TQMC-wh-n2 vs. wASTRID & abayes & 1X & deep & 200 & 33 / 10 / 7 & 5.899268e-05 & *** & MC"
#[1] "TQMC-wh-n2 vs. wASTRID & abayes & 1X & deep & 1000 & 30 / 6 / 14 & 5.731192e-04 & ** & MC"
#[1] "TQMC-wh-n2 vs. wASTRID & abayes & 1X & shallow & 50 & 30 / 13 / 7 & 1.018385e-03 & ** & "
#[1] "TQMC-wh-n2 vs. wASTRID & abayes & 1X & shallow & 200 & 26 / 13 / 11 & 1.804032e-02 & * & "
#[1] "TQMC-wh-n2 vs. wASTRID & abayes & 1X & shallow & 1000 & 26 / 12 / 12 & 9.458135e-03 & * & "
#[1] "TQMC-wh-n2 vs. wASTRID & abayes & 5X & deep & 50 & 31 / 7 / 12 & 2.961894e-05 & **** & MC"
#[1] "TQMC-wh-n2 vs. wASTRID & abayes & 5X & deep & 200 & 38 / 5 / 7 & 3.497277e-08 & ***** & MC"
#[1] "TQMC-wh-n2 vs. wASTRID & abayes & 5X & deep & 1000 & 43 / 2 / 5 & 5.250217e-09 & ***** & MC"
#[1] "TQMC-wh-n2 vs. wASTRID & abayes & 5X & shallow & 50 & 22 / 13 / 15 & 2.99532e-01 &  & "
#[1] "TQMC-wh-n2 vs. wASTRID & abayes & 5X & shallow & 200 & 12 / 20 / 18 & 9.461943e-01 &  & "
#[1] "TQMC-wh-n2 vs. wASTRID & abayes & 5X & shallow & 1000 & 14 / 7 / 29 & 2.197961e-01 &  & "
