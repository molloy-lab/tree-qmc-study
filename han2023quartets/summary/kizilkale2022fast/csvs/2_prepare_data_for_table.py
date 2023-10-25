import pandas
import math
import numpy
import sys


df = pandas.read_csv("data-all-error-and-timings_exfig4-6.csv")
df.fillna("NA", inplace=True)

cols = ["NCELL", "NMUT", "REPL",
        "TQMCn2xNQS", 
        "FASTRALxNQS",
        "SCISTREExNQS",
        "FASTMExNQS"]
rows = []

mthds = ["scistree_v1.2.0.6",
         "fastral",
         "fastme_v2.1.5",
         "treeqmcbip_v1.0.0_n2"]

ncxnms = [["n1000", "m300"], ["n300", "m300"], ["n300", "m1000"]]

for ncxnm in ncxnms:
    ncell, nmut = ncxnm

    repls = df[(df["NCELL"] == ncell) &
               (df["NMUT"] == nmut) &
               (df["MTHD"] == mthds[0])].REPL.values
    
    for repl in repls:
        print("%s %s %s" % (ncell, nmut, repl))
        xdf = df[(df["NCELL"] == ncell) &
                 (df["NMUT"] == nmut) &
                 (df["REPL"] == repl)]

        scistree = xdf[(xdf["MTHD"] == "scistree_v1.2.0.6")]
        fastral = xdf[(xdf["MTHD"] == "fastral")]
        fastme = xdf[(xdf["MTHD"] == "fastme_v2.1.5")]
        treeqmcbip = xdf[(xdf["MTHD"] == "treeqmcbip_v1.0.0_n2")]

        row = {}
        row["NCELL"] = ncell
        row["NMUT"] = nmut
        row["REPL"] = repl
        row["SCISTREExNQS"] = scistree.NQS.values[0]
        row["FASTRALxNQS"] = fastral.NQS.values[0]
        row["FASTMExNQS"] = fastme.NQS.values[0]
        row["TQMCn2xNQS"] = treeqmcbip.NQS.values[0]
        rows.append(row)

df = pandas.DataFrame(rows, columns=cols)

df.to_csv("data-for-table.csv",
          sep=',', na_rep="NA",header=True, index=False)
