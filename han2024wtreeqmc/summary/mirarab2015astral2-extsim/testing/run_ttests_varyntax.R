#!/usr/bin/env Rscript

# MANY NORMALITY FAILURES SO USE WILCOXON INSTEAD!!!

data <- read.csv("../csvs/data-varyntax-for-testing.csv")
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

ntests <- length(supps) * length(ngens) * 11 * length(mthds)
threshold0 <- 0.05
threshold_bonferoni <- threshold0 / ntests

writeme <- paste("alpha =", as.character(threshold0))
print(writeme)
print(paste("bonferroni alpha =", as.character(threshold_bonferoni), 
            "=", as.character(threshold0), "/", as.character(ntests)))

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
        tt <- t.test(mthd1, y = mthd2,
                     alternative="two.sided",
                     mu = 0, paired=TRUE, var.equal=FALSE,
                     conf.level=0.95)

        # Report
        if (ntreeqmc + nother == 0) {
            # all ties!
            tt$p.value <- "NA"
            stars = ""
            mc <- ""
        } else {
            # To be safe, check for normality
            #ggqqplot(diff)
            #hist(diff)
            st <- shapiro.test(diff)

            stars = ""
            if (tt$p.value < 0.000005) {
                stars = "*****"
            } else if (tt$p.value < 0.00005) {
                stars = "****"
            } else if (tt$p.value < 0.0005) {
                stars = "***"
            } else if (tt$p.value < 0.005) {
                stars = "**"
            } else if (tt$p.value < threshold0) {
                stars = "*"
            }
            mc <- ""
            if (tt$p.value < threshold_bonferoni) {
                mc <- "MC"
            }
            note <- ""
            if ((mean(diff) < 0) & (tt$p.value < threshold0)) {
                note <- paste("(", mthd, "better)")
            }
            if (st$p.value < 0.05) {
                note <- paste(note, "<- Failed normality test!")
            }
        }

        writeme <- paste("TQMC-wh-n2 vs.", mthd, "&",
                         supp, "&",
                         as.character(ntax), "&",
                         as.character(ngen), "&",
                         as.character(ntreeqmc), "/",
                         as.character(nother), "/",
                         as.character(ntie), "&",
                         format(tt$p.value, scientific = TRUE),
                         "&", stars, "&", mc, note)
        print(writeme)

            }
        }
    }
}

print("done")
exit()
