import numpy
import pandas
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import rc
import sys

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{helvet} \usepackage{sfmath}')

letters = [r"\textbf{a)}", 
           r"\textbf{b)}", 
           r"\textbf{c)}", 
           r"\textbf{d)}", 
           r"\textbf{e)}", 
           r"\textbf{f)}"]

# Tableau 20 colors in RGB.
darkblue = [(31, 119, 180), (174, 199, 232)]
orange = [(255, 127, 14), (255, 187, 120)]
green = [(44, 160, 44), (152, 223, 138)]
red = [(214, 39, 40), (255, 152, 150)]
purple = [(148, 103, 189), (197, 176, 213)]
brown = [(140, 86, 75), (196, 156, 148)]
pink = [(227, 119, 194), (247, 182, 210)]
gray = [(127, 127, 127), (199, 199, 199)]
yellow = [(188, 189, 34), (219, 219, 141)]
lightblue = [(23, 190, 207), (158, 218, 229)]

# Map RGB values to the [0, 1]
def map_rgb_to_01(tableau20):
    for i in range(len(tableau20)):
        r, g, b = tableau20[i]
        tableau20[i] = (r / 255., g / 255., b / 255.)

tableau20 = []
tableau20 += gray
tableau20 += brown
tableau20 += lightblue
tableau20 += darkblue
tableau20 += pink
tableau20 += purple
tableau20 += orange

map_rgb_to_01(tableau20)

z1 = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
z2 = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

fig = plt.figure(figsize=(10, 2))  # x, y
gs = gridspec.GridSpec(1, 5)  # nrows, ncols

ax00 = plt.subplot(gs[0, 0])  # Morel et al (2023)
ax01 = plt.subplot(gs[0, 1])  # vs. TREE-QMC shared (FNR)
ax02 = plt.subplot(gs[0, 2])  # vs. ASTRID          (FNR)
ax03 = plt.subplot(gs[0, 3])  # vs. ASTER           (FNR)
ax04 = plt.subplot(gs[0, 4])  # vs. ASTER           (LPP)



axs = [ax00, ax01, ax02, ax03, ax04]


# Read Morel et al (2023) CSVs, de-duplicate, and concatenate
print("Morel et al (2023)")
df1 = pandas.read_csv("morel2022asteroid-ilssim/csvs/data-varyblsc-for-testing.csv")
df2 = pandas.read_csv("morel2022asteroid-ilssim/csvs/data-varyngen-for-testing.csv")
df3 = pandas.read_csv("morel2022asteroid-ilssim/csvs/data-varymiss-for-testing.csv")
df4 = pandas.read_csv("morel2022asteroid-ilssim/csvs/data-varyntax-for-testing.csv")
df5 = pandas.read_csv("morel2022asteroid-ilssim/csvs/data-varynbps-for-testing.csv")
df6 = pandas.read_csv("morel2022asteroid-ilssim/csvs/data-varypsiz-for-testing.csv")

dfs = [df1]
for df in [df2, df3, df4, df5, df6]:
    df = df.drop(df[(df["NTAX"] == 50) &
                (df["NGEN"] == 1000) &
                (df["NBPS"] == 100) &
                (df["BLSC"] == 1.0) &
                (df["PSIZ"] == 50000000) &
                (df["MISS"] == 0.6)].index)
    dfs.append(df)

df = pandas.concat(dfs)

ax = ax00
ax.text(0.0, 0.5, "Morel et al. (2023)", fontsize=12)
ax.spines["bottom"].set_visible(False)
ax.spines["left"].set_visible(False)
ax.set_xticks([])
ax.set_yticks([])

# a) TREE-QMC (n2) vs. shared
ax = ax01
ax.set_title(letters[0], loc="left", fontsize=11)
x = df.TQMCn2refinedxFNR.values
y = df.TQMCn2sharedrefinedxFNR.values
ax.plot(x, y, '.', color=tableau20[1])
ax.plot(z1, z1, '-', color='r')
ax.set_aspect('equal', adjustable='box')
ax.set_xlabel("Norm RF")
ax.set_ylabel("TQMC (n2 shared)")
ax.set_xticks([0.0,  0.5, 1.0])
ax.set_yticks([0.0,  0.5, 1.0])


# b) TREE-QMC (n2) vs. ASTRID
ax = ax02
ax.set_title(letters[1], loc="left", fontsize=11)
y = df.WASTRIDxFNR.values
ax.plot(x, y, '.', color=tableau20[1])
ax.plot(z1, z1, '-', color='r')
ax.set_aspect('equal', adjustable='box')
ax.set_xlabel("Norm RF")
ax.set_ylabel("ASTRID")
ax.set_xticks([0.0,  0.5, 1.0])
ax.set_yticks([0.0,  0.5, 1.0])


# c) TREE-QMC (n2) vs. ASTRAL/ASTER
# SPECIES TREE ERROR - RF
ax = ax03
ax.set_title(letters[2], loc="left", fontsize=11)
y = df.ASTERxFNR.values
ax.plot(x, y, '.', color=tableau20[1])
ax.plot(z1, z1, '-', color='r')
ax.set_aspect('equal', adjustable='box')
ax.set_xlabel("Norm RF")
ax.set_ylabel("ASTRAL/ASTER")
ax.set_xticks([0.0,  0.5, 1.0])
ax.set_yticks([0.0,  0.5, 1.0])

diff = y - x                         # ASTER - TREE-QMC
n1 = len(numpy.where(diff == 0)[0])
n2 = len(numpy.where(diff > 0)[0])   # positive means TREE-QMC is better
n3 = len(numpy.where(diff < 0)[0])   # negative means ASTER is better
print("Norm RF --> TIE: %d, TQMC better: %d, ASTER better: %d, TOTAL: %d" % (n1, n2, n3, n1+n2+n3))


# QUARTET SCORE
# Didn't plot because not very illuminating (looks like straight line despite ASTRAL/ASTER being much better)
x = df.TQMCn2refinedxQS.values
y = df.ASTERxQS.values
diff = y - x                         # ASTER - TREE-QMC
n1 = len(numpy.where(diff == 0)[0])
n2 = len(numpy.where(diff < 0)[0])   # negative means TREE-QMC is better
n3 = len(numpy.where(diff > 0)[0])   # positive means ASTER is better
print("QS --> TIE: %d, TQMC better: %d, ASTER better: %d, TOTAL: %d" % (n1, n2, n3, n1+n2+n3))


# d) LOCAL PP
ax = ax04
ax.set_title(letters[3], loc="left", fontsize=11)
x = df.TQMCn2refinedxAVGLPP.values
y = df.ASTERxAVGLPP.values
ax.plot(x, y, '.', color=tableau20[1])
ax.plot(z2, z2, '-', color='r')
ax.set_aspect('equal', adjustable='box')
ax.set_xlabel("Local PP")
ax.set_ylabel("ASTRAL/ASTER")
ax.set_xticks([0.5, 0.75, 1.0])
ax.set_yticks([0.5, 0.75, 1.0])

diff = y - x                         # ASTER - TREE-QMC
n1 = len(numpy.where(diff == 0)[0])
n2 = len(numpy.where(diff < 0)[0])   # negative means TREE-QMC is better
n3 = len(numpy.where(diff > 0)[0])   # positive means ASTER is better
print("Mean Local PP --> TIE: %d, TQMC better: %d, ASTER better: %d, TOTAL: %d" % (n1, n2, n3, n1+n2+n3))


for ax in axs:
    ax.tick_params(axis='x', labelsize=9)
    ax.tick_params(axis='y', labelsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


gs.tight_layout(fig, rect=[0, 0, 1, 1])
plt.savefig("test.pdf", format='pdf', dpi=300)


# Morel et al (2023)
# Norm RF       --> TIE: 178, TQMC better: 801, ASTER better: 321,  TOTAL: 1300
# QS            --> TIE: 21,  TQMC better: 13,  ASTER better: 1266, TOTAL: 1300
# Mean Local PP --> TIE: 9,   TQMC better: 725, ASTER better: 566,  TOTAL: 1300









