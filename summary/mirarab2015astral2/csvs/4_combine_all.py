import pandas
import numpy
import sys

# Combine all CSVs
df1 = pandas.read_csv("data-lt200tax-error-and-timings.csv", 
                      keep_default_na=False)
df2 = pandas.read_csv("data-eq200tax-error-and-timings.csv", 
                      keep_default_na=False)
df3 = pandas.read_csv("data-gt200tax-error-and-timings.csv", 
                      keep_default_na=False)

new_df = pandas.concat([df1, df2, df3])
new_df.to_csv("data-all-error-and-timings.csv",
              sep=',', na_rep="NA", header=True, index=False)
