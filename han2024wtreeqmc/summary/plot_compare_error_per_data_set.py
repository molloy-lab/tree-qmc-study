import numpy
import pandas
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import rc
import sys

# TO-DO Add something that says, north of line equals TREE-QMC is better

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{helvet} \usepackage{sfmath}')

letters = [r"\textbf{a)}", 
           r"\textbf{b)}", 
           r"\textbf{c)}", 
           r"\textbf{d)}", 
           r"\textbf{e)}", 
           r"\textbf{f)}",
           r"\textbf{g)}",
           r"\textbf{h)}",
           r"\textbf{i)}",
           r"\textbf{j)}",
           r"\textbf{k)}",
           r"\textbf{l)}",
           r"\textbf{m)}",
           r"\textbf{n)}"]

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

fig = plt.figure(figsize=(10, 4))  # x, y
gs = gridspec.GridSpec(2, 5)  # nrows, ncols (reversed)

ax00 = plt.subplot(gs[0, 0])  # Morel et al (2023)
ax01 = plt.subplot(gs[0, 1])  # vs. TREE-QMC shared (FNR)
ax02 = plt.subplot(gs[0, 2])  # vs. ASTRID          (FNR)
ax03 = plt.subplot(gs[0, 3])  # vs. ASTER           (FNR)
ax04 = plt.subplot(gs[0, 4])  # vs. ASTER           (LPP)

ax10 = plt.subplot(gs[1, 0])  # Morel et al (2023)
ax11 = plt.subplot(gs[1, 1])  # vs. TREE-QMC shared (FNR)
ax12 = plt.subplot(gs[1, 2])  # vs. ASTRID          (FNR)
ax13 = plt.subplot(gs[1, 3])  # vs. ASTER           (FNR)
ax14 = plt.subplot(gs[1, 4])  # vs. ASTER           (LPP)


axs = [ax00, ax01, ax02, ax03, ax04,
       ax10, ax11, ax12, ax13, ax14]

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

# e, f, g)
x = df.TQMCn2refinedxFNR.values

# a) TREE-QMC (n2) vs. shared
ax = ax01
ax.set_title(letters[0], loc="left", fontsize=11)
y = df.TQMCn2sharedrefinedxFNR.values
ax.plot(x, y, '.', color=tableau20[1])
ax.plot([0, 1], [0, 1], '-', color='r')
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
ax.plot([0, 1], [0, 1], '-', color='r')
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
ax.plot([0, 1], [0, 1], '-', color='r')
ax.set_aspect('equal', adjustable='box')
ax.set_xlabel("Norm RF")
ax.set_ylabel("ASTRAL/ASTER")
ax.set_xticks([0.0,  0.5, 1.0])
ax.set_yticks([0.0,  0.5, 1.0])

diff = y - x                         # ASTER - TREE-QMC
n1 = len(numpy.where(diff == 0)[0])
n2 = len(numpy.where(diff > 0)[0])   # positive means TREE-QMC is better
n3 = len(numpy.where(diff < 0)[0])   # negative means ASTER is better
print("      Norm RF --> TIE: %d, TQMC better: %d, ASTER better: %d, TOTAL: %d" % (n1, n2, n3, n1+n2+n3))


# QUARTET SCORE
# Didn't plot because not very illuminating (looks like straight line despite ASTRAL/ASTER being much better)
x = df.TQMCn2refinedxQS.values
y = df.ASTERxQS.values
diff = y - x                         # ASTER - TREE-QMC
n1 = len(numpy.where(diff == 0)[0])
n2 = len(numpy.where(diff < 0)[0])   # negative means TREE-QMC is better
n3 = len(numpy.where(diff > 0)[0])   # positive means ASTER is better
print("           QS --> TIE: %d, TQMC better: %d, ASTER better: %d, TOTAL: %d" % (n1, n2, n3, n1+n2+n3))


# d) LOCAL PP
ax = ax04
ax.set_title(letters[3], loc="left", fontsize=11)
x = df.TQMCn2refinedxAVGLPP.values
y = df.ASTERxAVGLPP.values
ax.plot(x, y, '.', color=tableau20[1])
ax.plot([0.5, 1.0], [0.5, 1.0], '-', color='r')
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

# Morel et al (2023)
# Norm RF       --> TIE: 178, TQMC better: 801, ASTER better: 321,  TOTAL: 1300
# QS            --> TIE: 21,  TQMC better: 13,  ASTER better: 1266, TOTAL: 1300
# Mean Local PP --> TIE: 9,   TQMC better: 725, ASTER better: 566,  TOTAL: 1300

###########

# Read Zhang et al (2018) CSVs and remove bs support (focus on abayes)
print("\nZhang et al (2022)")
df = pandas.read_csv("zhang2022weighting-gteesim/csvs/data-for-testing.csv")
df = df.drop(df[(df["SUPP"] == "bs")].index)

ax = ax10
ax.text(0.0, 0.5, "Zhang et al. (2018)\n w/ abayes support", fontsize=12)
ax.spines["bottom"].set_visible(False)
ax.spines["left"].set_visible(False)
ax.set_xticks([])
ax.set_yticks([])

# e, f, g)
x = df.TQMCwhn2xSERF.values

# e) TREE-QMC (wh, n2) vs. unweighted TREE-QMC (wf, n2)
ax = ax11
ax.set_title(letters[4], loc="left", fontsize=11)
y = df.TQMCn2xSERF.values
ax.plot(x, y, '.', color=tableau20[1])
ax.plot([0.0, 0.5], [0.0, 0.5], '-', color='r')
ax.set_aspect('equal', adjustable='box')
ax.set_xlabel("Norm RF")
ax.set_ylabel("TQMC (unweighted)")
ax.set_xticks([0.0, 0.25, 0.5])
ax.set_yticks([0.0, 0.25, 0.5])

# f) TREE-QMC (wh, n2) vs. ASTRID (ws)
ax = ax12
ax.set_title(letters[5], loc="left", fontsize=11)
y = df.WASTRIDxSERF.values
ax.plot(x, y, '.', color=tableau20[1])
ax.plot([0.0, 0.4], [0.0, 0.4], '-', color='r')
ax.set_aspect('equal', adjustable='box')
ax.set_xlabel("Norm RF")
ax.set_ylabel("ASTRID-ws")
ax.set_xticks([0.0,  0.2, 0.4])
ax.set_yticks([0.0,  0.2, 0.4])

# g) TREE-QMC (n2) vs. ASTRAL/ASTER
# SPECIES TREE ERROR - RF
ax = ax13
ax.set_title(letters[6], loc="left", fontsize=11)
y = df.ASTERHxSERF.values
ax.plot(x, y, '.', color=tableau20[1])
ax.plot([0, 0.4], [0, 0.4], '-', color='r')
ax.set_aspect('equal', adjustable='box')
ax.set_xlabel("Norm RF")
ax.set_ylabel("ASTRAL/ASTER (wh)")
ax.set_xticks([0.0,  0.2, 0.4])
ax.set_yticks([0.0,  0.2, 0.4])

diff = y - x                         # ASTER - TREE-QMC
n1 = len(numpy.where(diff == 0)[0])
n2 = len(numpy.where(diff > 0)[0])   # positive means TREE-QMC is better
n3 = len(numpy.where(diff < 0)[0])   # negative means ASTER is better
print("      Norm RF --> TIE: %d, TQMC better: %d, ASTER better: %d, TOTAL: %d" % (n1, n2, n3, n1+n2+n3))


# QUARTET SCORE
# Didn't plot because not very illuminating (looks like straight line despite ASTRAL/ASTER being much better)
x = df.TQMCwhn2xQS.values
y = df.ASTERHxQS.values
diff = y - x                         # ASTER - TREE-QMC
n1 = len(numpy.where(diff == 0)[0])
n2 = len(numpy.where(diff < 0)[0])   # negative means TREE-QMC is better
n3 = len(numpy.where(diff > 0)[0])   # positive means ASTER is better
print("           QS --> TIE: %d, TQMC better: %d, ASTER better: %d, TOTAL: %d" % (n1, n2, n3, n1+n2+n3))


# h) LOCAL PP
ax = ax14
ax.set_title(letters[3], loc="left", fontsize=11)
x = df.TQMCwhn2xAVGLPP.values
y = df.ASTERHxAVGLPP.values
ax.plot(x, y, '.', color=tableau20[1])
ax.plot([0.8, 1.0], [0.8, 1.0], '-', color='r')
ax.set_aspect('equal', adjustable='box')
ax.set_xlabel("Local PP")
ax.set_ylabel("ASTRAL/ASTER (wh)")
ax.set_xticks([0.8, 0.9, 1.0])
ax.set_yticks([0.8, 0.9, 1.0])

diff = y - x                         # ASTER - TREE-QMC
n1 = len(numpy.where(diff == 0)[0])
n2 = len(numpy.where(diff < 0)[0])   # negative means TREE-QMC is better
n3 = len(numpy.where(diff > 0)[0])   # positive means ASTER is better
print("Mean Local PP --> TIE: %d, TQMC better: %d, ASTER better: %d, TOTAL: %d" % (n1, n2, n3, n1+n2+n3))


###########




###########

for ax in axs:
    ax.tick_params(axis='x', labelsize=9)
    ax.tick_params(axis='y', labelsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


gs.tight_layout(fig, rect=[0, 0, 1, 1])
plt.savefig("test.pdf", format='pdf', dpi=300)













