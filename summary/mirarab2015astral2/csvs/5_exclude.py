import pandas
import numpy
import sys


df = pandas.read_csv("data-all-error-and-timings.csv", keep_default_na=False)

ntaxs = [10, 50, 100, 200, 500, 1000]
gtres = ["estimatedgenetre", "truegenetrees"]
ngens = [250, 1000]
dfs = []

for ntax in ntaxs:
    if ntax == 200:
        strhts = [10000000, 2000000, 500000]
        srates = [0.0000001, 0.000001]
    else:
        strhts = [2000000]
        srates = [0.000001]
    for strht in strhts:
        for srate in srates:
            for ngen in ngens:
                for gtre in gtres:
                    repls = list(range(1, 51))

                    if gtre == "estimatedgenetre":
                        if ntax == 10:
                            repls.remove(41)
                        elif ntax == 50:
                            repls.remove(21)
                            repls.remove(41)
                        elif ntax == 100:
                            repls.remove(8)
                            repls.remove(47)
                        elif ntax == 200:
                            if (strht == 500000) and (srate == 0.000001):
                                repls.remove(8)
                                repls.remove(15)
                                repls.remove(49)
                        elif ntax == 1000:
                            if ngen == 1000:
                                # EKM manually confirmed ASTRAL-III didn't 
                                # finish in 18 hours on these three replicate 
                                # data sets for one model condition.
                                repls.remove(6)
                                repls.remove(8)
                                repls.remove(38)

                    for repl in repls:
                        xdf = df[(df["NTAX"] == ntax) &
                                 (df["STRHT"] == strht) &
                                 (df["SRATE"] == srate) &
                                 (df["NGEN"] == ngen) &
                                 (df["GTRE"] == gtre) &  
                                 (df["REPL"] == repl)]
                        dfs.append(xdf)
    
new_df = pandas.concat(dfs)

new_df.to_csv("data-all-error-and-timings-ex.csv",
              sep=',', na_rep="NA", header=True, index=False)
