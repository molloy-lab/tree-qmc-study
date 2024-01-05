#!/usr/bin/env Rscript

data <- read.csv("../csvs/data-for-testing-varyntax.csv")
data$NTAX <- as.factor(data$NTAX)
data$HGHT <- as.factor(data$HGHT)
data$RATE <- as.factor(data$RATE)
data$SUPP <- as.factor(data$SUPP)
data$NGEN <- as.factor(data$NGEN)
data$REPL <- as.factor(data$REPL)

#data <- subset(data, HGHT == "2000000")
#data <- subset(data, RATE == "1e-06")

supps <- c("abayes")
ngens <- c(50, 200, 1000)

ntaxs <- c(10, 50, 100, 200, 500, 1000)

mthds <- c("CAML", "ASTER-h", "wASTRID")

ntest <- length(supps) * length(ngens) * 
         length(ntaxs) * length(mthds)

threshold_bonferoni <- 0.05 / ntest

for (mthd in mthds) {
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

		#if (0 %in% fndf$DIFF) {
		#	print("Found a tie!")
		#}

        # Count ties
        ntreeqmc <- sum(fndf$DIFF > 0)
		nother <- sum(fndf$DIFF < 0)
        ntie <- sum(fndf$DIFF == 0)


		
        # Compute average difference
        #avgdiffwtie <- mean(ydf$DIFF)
        #vardiffwtie <- sd(ydf$DIFF)
        #meddiffwtie <- median(ydf$DIFF)

        #avgdiff <- mean(zdf$DIFF)
        #meddiff <- median(zdf$DIFF)

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
                         as.character(ntax), "&",
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


#[1] "TQMC-wh-n2 vs. CAML & abayes & 10 & 50 & 8 / 6 / 35 & 4.725395e-01 &  & "
#[1] "TQMC-wh-n2 vs. CAML & abayes & 10 & 200 & 4 / 2 / 43 & 3.223252e-01 &  & "
#[1] "TQMC-wh-n2 vs. CAML & abayes & 10 & 1000 & 2 / 1 / 45 & 4.199676e-01 &  & "
#[1] "TQMC-wh-n2 vs. CAML & abayes & 50 & 50 & 28 / 9 / 10 & 1.166386e-04 & *** & MC"
#[1] "TQMC-wh-n2 vs. CAML & abayes & 50 & 200 & 29 / 12 / 6 & 6.405316e-04 & ** & MC"
#[1] "TQMC-wh-n2 vs. CAML & abayes & 50 & 1000 & 21 / 3 / 23 & 1.761287e-03 & ** & "
#[1] "TQMC-wh-n2 vs. CAML & abayes & 100 & 50 & 38 / 5 / 5 & 3.302366e-08 & ***** & MC"
#[1] "TQMC-wh-n2 vs. CAML & abayes & 100 & 200 & 32 / 9 / 7 & 1.917265e-03 & ** & "
#[1] "TQMC-wh-n2 vs. CAML & abayes & 100 & 1000 & 23 / 13 / 12 & 4.272706e-02 & * & "
#[1] "TQMC-wh-n2 vs. CAML & abayes & 200 & 50 & 44 / 3 / 3 & 5.084971e-12 & ***** & MC"
#[1] "TQMC-wh-n2 vs. CAML & abayes & 200 & 200 & 40 / 4 / 6 & 3.090232e-09 & ***** & MC"
#[1] "TQMC-wh-n2 vs. CAML & abayes & 200 & 1000 & 33 / 8 / 9 & 1.054732e-03 & ** & "
#[1] "TQMC-wh-n2 vs. CAML & abayes & 500 & 50 & 46 / 3 / 1 & 8.033952e-16 & ***** & MC"
#[1] "TQMC-wh-n2 vs. CAML & abayes & 500 & 200 & 42 / 7 / 1 & 1.785721e-09 & ***** & MC"
#[1] "TQMC-wh-n2 vs. CAML & abayes & 500 & 1000 & 36 / 11 / 3 & 3.079315e-02 & * & "
#[1] "TQMC-wh-n2 vs. CAML & abayes & 1000 & 50 & 45 / 5 / 0 & 3.463452e-09 & ***** & MC"
#[1] "TQMC-wh-n2 vs. CAML & abayes & 1000 & 200 & 42 / 8 / 0 & 4.544811e-04 & *** & MC"

#[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 10 & 50 & 1 / 1 / 47 & 1e+00 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 10 & 200 & 0 / 1 / 48 & 3.223252e-01 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 10 & 1000 & 0 / 0 / 48 & NA &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 50 & 50 & 11 / 8 / 28 & 3.22542e-01 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 50 & 200 & 8 / 3 / 36 & 1.331484e-01 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 50 & 1000 & 5 / 3 / 39 & 4.854166e-01 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 100 & 50 & 17 / 8 / 23 & 2.4382e-02 & * & "
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 100 & 200 & 15 / 7 / 26 & 8.637687e-02 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 100 & 1000 & 8 / 10 / 30 & 6.521362e-01 &  & "
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 200 & 50 & 22 / 17 / 11 & 1.697341e-02 & * & "
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 200 & 200 & 25 / 10 / 15 & 1.244558e-02 & * & "
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 200 & 1000 & 17 / 5 / 28 & 6.8334e-03 & * & "
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 500 & 50 & 38 / 10 / 2 & 1.014811e-07 & ***** & MC"
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 500 & 200 & 40 / 6 / 4 & 1.113474e-08 & ***** & MC"
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 500 & 1000 & 36 / 8 / 6 & 2.573416e-05 & **** & MC"
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 1000 & 50 & 44 / 2 / 4 & 1.201765e-06 & ***** & MC"
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 1000 & 200 & 38 / 4 / 8 & 7.387226e-06 & **** & MC"
#[1] "TQMC-wh-n2 vs. ASTER-h & abayes & 1000 & 1000 & 44 / 2 / 4 & 1.071559e-04 & *** & MC"

#[1] "TQMC-wh-n2 vs. wASTRID & abayes & 10 & 50 & 1 / 2 / 46 & 5.690627e-01 &  & "
#[1] "TQMC-wh-n2 vs. wASTRID & abayes & 10 & 200 & 0 / 1 / 48 & 3.223252e-01 &  & "
#[1] "TQMC-wh-n2 vs. wASTRID & abayes & 10 & 1000 & 0 / 0 / 48 & NA &  & "
#[1] "TQMC-wh-n2 vs. wASTRID & abayes & 50 & 50 & 12 / 13 / 22 & 7.848466e-01 &  & "
#[1] "TQMC-wh-n2 vs. wASTRID & abayes & 50 & 200 & 18 / 3 / 26 & 2.530514e-03 & ** & "
#[1] "TQMC-wh-n2 vs. wASTRID & abayes & 50 & 1000 & 6 / 3 / 38 & 3.22542e-01 &  & "
#[1] "TQMC-wh-n2 vs. wASTRID & abayes & 100 & 50 & 21 / 14 / 13 & 1.604197e-01 &  & "
#[1] "TQMC-wh-n2 vs. wASTRID & abayes & 100 & 200 & 19 / 13 / 16 & 9.167777e-02 &  & "
#[1] "TQMC-wh-n2 vs. wASTRID & abayes & 100 & 1000 & 14 / 7 / 27 & 2.750208e-02 & * & "
#[1] "TQMC-wh-n2 vs. wASTRID & abayes & 200 & 50 & 30 / 13 / 7 & 1.018385e-03 & ** & "
#[1] "TQMC-wh-n2 vs. wASTRID & abayes & 200 & 200 & 26 / 13 / 11 & 1.804032e-02 & * & "
#[1] "TQMC-wh-n2 vs. wASTRID & abayes & 200 & 1000 & 26 / 12 / 12 & 9.458135e-03 & * & "
#[1] "TQMC-wh-n2 vs. wASTRID & abayes & 500 & 50 & 35 / 13 / 2 & 2.053157e-04 & *** & MC"
#[1] "TQMC-wh-n2 vs. wASTRID & abayes & 500 & 200 & 41 / 7 / 2 & 5.374896e-07 & ***** & MC"
#[1] "TQMC-wh-n2 vs. wASTRID & abayes & 500 & 1000 & 29 / 9 / 12 & 1.040949e-04 & *** & MC"
#[1] "TQMC-wh-n2 vs. wASTRID & abayes & 1000 & 50 & 43 / 6 / 1 & 1.877626e-06 & ***** & MC"
#[1] "TQMC-wh-n2 vs. wASTRID & abayes & 1000 & 200 & 46 / 3 / 1 & 6.873988e-06 & **** & MC"
#[1] "TQMC-wh-n2 vs. wASTRID & abayes & 1000 & 1000 & 41 / 6 / 3 & 2.004529e-04 & *** & MC"
