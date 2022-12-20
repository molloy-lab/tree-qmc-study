import pandas
import math
import numpy
import sys


df = pandas.read_csv("data-all-error-and-timings-ex.csv")
df.fillna("NA", inplace=True)

cols = ["NTAX", "STRHT", "SRATE", "REPL",
        "GTRE", "NGEN", "NODE",
        "TQMCn2xSEFN",  "TQMCn2xSERF",  "TQMCn2xNQSC",  "TQMCn2xSECS",
        "FASTRALxSEFN", "FASTRALxSERF", "FASTRALxNQSC", "FASTRALxSECS",
        "ASTRAL3xSEFN", "ASTRAL3xSERF", "ASTRAL3xNQSC", "ASTRAL3xSECS"]
rows = []

gtre = "estimatedgenetre"
ntaxs = [10, 50, 100, 200, 500, 1000]
ngens = [250, 1000]
mthds = ["treeqmc_n0_v1.0.0",
         "treeqmc_n1_v1.0.0",
         "treeqmc_n2_v1.0.0",
         
         "fastral"]
mthds = ["treeqmc_n2_v1.0.0",
         "fastral",
         "astral_3_v5.7.7",]

for ntax in ntaxs:
    repls = list(range(1, 51))
    if ntax == 200:
        strhts = [10000000, 2000000, 500000]
        srates = [0.0000001, 0.000001]
    else:
        strhts = [2000000]
        srates = [0.000001]
    for strht in strhts:
        for srate in srates:
            for ngen in ngens:
                repls = df[(df["NTAX"] == ntax) &
                           (df["STRHT"] == strht) &
                           (df["SRATE"] == srate) &
                           (df["GTRE"] == gtre) &
                           (df["NGEN"] == ngen) &
                           (df["MTHD"] == mthds[0])].REPL.values
                for repl in repls:
                    print("%d %d %g %d %d" % (ntax, strht, srate, repl, ngen))
                    xdf = df[(df["NTAX"] == ntax) &
                             (df["STRHT"] == strht) &
                             (df["SRATE"] == srate) &
                             (df["REPL"] == repl) & 
                             (df["GTRE"] == gtre) &
                             (df["NGEN"] == ngen)]

                    #tqmc0 = xdf[(xdf["MTHD"] == "treeqmc_n0_v1.0.0")]
                    #tqmc1 = xdf[(xdf["MTHD"] == "treeqmc_n1_v1.0.0")]
                    tqmc2 = xdf[(xdf["MTHD"] == "treeqmc_n2_v1.0.0")]
                    fastral = xdf[(xdf["MTHD"] == "fastral")]
                    astral = xdf[(xdf["MTHD"] == "astral_3_v5.7.7")]

                    node1 = tqmc2.NODE.values[0]
                    for mthd in mthds:
                        node2 = xdf[(xdf["MTHD"] == mthd)].NODE.values[0]
                        if (node1 != node2) and (node2 != "NA"):
                            print(xdf)
                            sys.exit("Node error!")

                    row = {}
                    row["NTAX"] = ntax
                    row["STRHT"] = strht
                    row["SRATE"] = srate
                    row["REPL"] = repl
                    row["GTRE"] = gtre
                    row["NGEN"] = ngen
                    row["NODE"] = node1

                    #row["TQMCn0xSEFN"] = tqmc0.SEFN.values[0]
                    #row["TQMCn0xSERF"] = tqmc0.SERF.values[0]
                    #row["TQMCn0xNQSC"] = tqmc0.NQSC.values[0]
                    #row["TQMCn0xSECS"] = tqmc0.SECS.values[0]

                    #row["TQMCn1xSEFN"] = tqmc1.SEFN.values[0]
                    #row["TQMCn1xSERF"] = tqmc1.SERF.values[0]
                    #row["TQMCn1xNQSC"] = tqmc1.NQSC.values[0]
                    #row["TQMCn1xSECS"] = tqmc1.SECS.values[0]

                    row["TQMCn2xSEFN"] = tqmc2.SEFN.values[0]
                    row["TQMCn2xSERF"] = tqmc2.SERF.values[0]
                    row["TQMCn2xNQSC"] = tqmc2.NQSC.values[0]
                    row["TQMCn2xSECS"] = tqmc2.SECS.values[0]

                    row["FASTRALxSEFN"] = fastral.SEFN.values[0]
                    row["FASTRALxSERF"] = fastral.SERF.values[0]
                    row["FASTRALxNQSC"] = fastral.NQSC.values[0]
                    row["FASTRALxSECS"] = fastral.SECS.values[0]

                    row["ASTRAL3xSEFN"] = astral.SEFN.values[0]
                    row["ASTRAL3xSERF"] = astral.SERF.values[0]
                    row["ASTRAL3xNQSC"] = astral.NQSC.values[0]
                    row["ASTRAL3xSECS"] = astral.SECS.values[0]

                    rows.append(row)

df = pandas.DataFrame(rows, columns=cols)

#xdf = df[(df["NTAX"] == 1000) & (df["NGEN"] == 1000)]
#print(xdf.FASTRALxSEFN.values)
#print(xdf.TQMCn2xSEFN.values)

df.to_csv("data-for-testing.csv",
          sep=',', na_rep="NA",header=True, index=False)
