import numpy
import pandas
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import rc
from matplotlib import cm
from matplotlib.colors import Normalize
from scipy.interpolate import interpn
import sys

def density_scatter(x, y, ax, sort=True, bins=20, **kwargs):
    """
    Scatter plot colored by 2d histogram
    https://stackoverflow.com/questions/20105364/how-can-i-make-a-scatter-plot-colored-by-density
    """
    data , x_e, y_e = numpy.histogram2d(x, y, bins=bins, density=True)
    z = interpn( ( 0.5*(x_e[1:] + x_e[:-1]) , 0.5*(y_e[1:] + y_e[:-1]) ), 
                 data , 
                 numpy.vstack([x,y]).T, 
                 method = "splinef2d", 
                 bounds_error = False)

    #To be sure to plot all data
    z[numpy.where(numpy.isnan(z))] = 0.0

    # Sort the points by density, so that the densest points are plotted last
    if sort :
        idx = z.argsort()
        x, y, z = x[idx], y[idx], z[idx]

    ax.scatter( x, y, c=z, marker='.', **kwargs )

    #norm = Normalize(vmin = numpy.min(z), vmax = numpy.max(z))
    #cbar = fig.colorbar(cm.ScalarMappable(norm = norm), ax=ax)
    #cbar.ax.set_ylabel('Density')

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

fig = plt.figure(figsize=(10, 6))  # x, y
gs = gridspec.GridSpec(3, 5)       # nrows, ncols (reversed)

ax00 = plt.subplot(gs[0, 0])
ax01 = plt.subplot(gs[0, 1])
ax02 = plt.subplot(gs[0, 2])
ax03 = plt.subplot(gs[0, 3])
ax04 = plt.subplot(gs[0, 4])

ax10 = plt.subplot(gs[1, 0])
ax11 = plt.subplot(gs[1, 1])
ax12 = plt.subplot(gs[1, 2])
ax13 = plt.subplot(gs[1, 3])
ax14 = plt.subplot(gs[1, 4])

ax20 = plt.subplot(gs[2, 0])
ax21 = plt.subplot(gs[2, 1])
ax22 = plt.subplot(gs[2, 2])
ax23 = plt.subplot(gs[2, 3])
ax24 = plt.subplot(gs[2, 4])

axs = [ax00, ax01, ax02, ax03, ax04,
       ax10, ax11, ax12, ax13, ax14,
       ax20, ax21, ax22, ax23, ax24]

########### Read Morel et al (2023)

# Read CSVs, de-duplicate, and concatenate
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
ax.text(0.0, 0.5, "Data sets from\nMorel et al. (2023)", fontsize=12)
ax.spines["bottom"].set_visible(False)
ax.spines["left"].set_visible(False)
ax.set_xticks([])
ax.set_yticks([])

# e, f, g)
x = df.TQMCn2refinedxFNR.values * 100

# a) TREE-QMC (n2) vs. shared
name = "shared"
y = df.TQMCn2sharedrefinedxFNR.values * 100

ax = ax01
ax.set_title(letters[0], loc="left", fontsize=11)
density_scatter(x, y, ax) 
#ax.plot(x, y, '.', color=tableau20[1])
ax.plot([0, 100], [0, 100], '-', color='r')
ax.set_aspect('equal', adjustable='box')
ax.set_ylabel("(shared n2)")
ax.set_xlabel("TREE-QMC (n2)")
ax.set_xticks([0,  50, 100])
ax.set_yticks([0,  50, 100])


diff = y - x
n1 = len(numpy.where(diff == 0)[0])
n2 = len(numpy.where(diff > 0)[0])   # negative means TREE-QMC is better
n3 = len(numpy.where(diff < 0)[0])   # positive means other method is better
print("      Norm RF --> TIE: %d,  TQMC better: %d, %s better: %d, TOTAL: %d" % (n1, n2, name, n3, n1+n2+n3))

# b) TREE-QMC (n2) vs. ASTRID
name = "ASTRID"
y = df.WASTRIDxFNR.values * 100

ax = ax02
ax.set_title(letters[1], loc="left", fontsize=11)
density_scatter(x, y, ax) 
#ax.plot(x, y, '.', color=tableau20[1])
ax.plot([0, 100], [0, 100], '-', color='r')
ax.set_aspect('equal', adjustable='box')
ax.set_ylabel(name)
ax.set_xlabel("TREE-QMC (n2)")
ax.set_xticks([0,  50, 100])
ax.set_yticks([0,  50, 100])

diff = y - x
n1 = len(numpy.where(diff == 0)[0])
n2 = len(numpy.where(diff > 0)[0])   # positive means TREE-QMC is better
n3 = len(numpy.where(diff < 0)[0])   # negative means other method is better
print("      Norm RF --> TIE: %d, TQMC better: %d, %s better: %d, TOTAL: %d" % (n1, n2, name, n3, n1+n2+n3))


name = "ASTEROID"
y = df.ASTEROIDxFNR.values * 100
diff = y - x
n1 = len(numpy.where(diff == 0)[0])
n2 = len(numpy.where(diff > 0)[0])   # positive means TREE-QMC is better
n3 = len(numpy.where(diff < 0)[0])   # negative means other method is better
print("      Norm RF --> TIE: %d, TQMC better: %d, %s better: %d, TOTAL: %d" % (n1, n2, name, n3, n1+n2+n3))


# c) TREE-QMC (n2) vs. ASTRAL/ASTER - SPECIES TREE ERROR - RF
name = "ASTRAL"
y = df.ASTERxFNR.values * 100

ax = ax03
ax.set_title(letters[2], loc="left", fontsize=11)
density_scatter(x, y, ax)
#ax.plot(x, y, '.', color=tableau20[1])
ax.plot([0, 100], [0, 100], '-', color='r')
ax.set_aspect('equal', adjustable='box')
ax.set_ylabel(name)
ax.set_xlabel("TREE-QMC (n2)")
ax.set_xticks([0,  50, 100])
ax.set_yticks([0,  50, 100])

diff = y - x
n1 = len(numpy.where(diff == 0)[0])
n2 = len(numpy.where(diff > 0)[0])   # positive means TREE-QMC is better
n3 = len(numpy.where(diff < 0)[0])   # negative means other method is better
print("      Norm RF --> TIE: %d, TQMC better: %d, %s better: %d, TOTAL: %d" % (n1, n2, name, n3, n1+n2+n3))


# QUARTET SCORE
# Looks like straight line despite ASTRAL/ASTER being much better
x = df.TQMCn2refinedxQS.values
y = df.ASTERxQS.values

diff = y - x
n1 = len(numpy.where(diff == 0)[0])
n2 = len(numpy.where(diff < 0)[0])   # negative means TREE-QMC is better
n3 = len(numpy.where(diff > 0)[0])   # positive means other method is better
print("           QS --> TIE: %d,  TQMC better: %d, %s better: %d, TOTAL: %d" % (n1, n2, name, n3, n1+n2+n3))


# d) TREE-QMC (n2) vs. ASTRAL/ASTER - LOCAL PP
x = df.TQMCn2refinedxAVGLPP.values
y = df.ASTERxAVGLPP.values

ax = ax04
ax.set_title(letters[3], loc="left", fontsize=11)
density_scatter(x, y, ax)
#ax.plot(x, y, '.', color=tableau20[1])
ax.plot([0.5, 1.0], [0.5, 1.0], '-', color='r')
ax.set_aspect('equal', adjustable='box')
ax.set_ylabel(name)
ax.set_xlabel("TREE-QMC (n2)")
ax.set_xticks([0.5, 0.75, 1.0])
ax.set_yticks([0.5, 0.75, 1.0])

diff = y - x
n1 = len(numpy.where(diff == 0)[0])
n2 = len(numpy.where(diff < 0)[0])   # negative means TREE-QMC is better
n3 = len(numpy.where(diff > 0)[0])   # positive means other method is better
print("Mean Local PP --> TIE: %d,   TQMC better: %d, %s better: %d, TOTAL: %d" % (n1, n2, name, n3, n1+n2+n3))

############ Zhang et al (2018)

# Read CSVs and remove bs support (focus on abayes)
print("\nZhang et al (2022)")
df = pandas.read_csv("zhang2022weighting-gteesim/csvs/data-for-testing.csv")
df = df.drop(df[(df["SUPP"] == "bs")].index)

ax = ax10
ax.text(0.0, 0.5, "Data sets from\nZhang et al. (2018)\nw/ abayes support", fontsize=12)
ax.spines["bottom"].set_visible(False)
ax.spines["left"].set_visible(False)
ax.set_xticks([])
ax.set_yticks([])

# e, f, g)
x = df.TQMCwhn2xSERF.values * 100

# e) TREE-QMC-wh (n2) vs. unweighted TREE-QMC (n2)
name = "unweighted (n2)"
y = df.TQMCn2xSERF.values * 100

ax = ax11
ax.set_title(letters[4], loc="left", fontsize=11)
density_scatter(x, y, ax)
#ax.plot(x, y, '.', color=tableau20[1])
ax.plot([0, 40], [0, 40], '-', color='r')
ax.set_aspect('equal', adjustable='box')
ax.set_ylabel(name)
ax.set_xlabel("TREE-QMC-wh (n2)")
ax.set_xticks([0, 20, 40])
ax.set_yticks([0, 20, 40])

diff = y - x
n1 = len(numpy.where(diff == 0)[0])
n2 = len(numpy.where(diff > 0)[0])   # positive means TREE-QMC is better
n3 = len(numpy.where(diff < 0)[0])   # negative means other method is better
print("      Norm RF --> TIE: %d, TQMC better: %d, %s better: %d, TOTAL: %d" % (n1, n2, name, n3, n1+n2+n3))

# f) TREE-QMC-wh (n2) vs. ASTRID-ws
name = "ASTRID-ws"
y = df.WASTRIDxSERF.values * 100

ax = ax12
ax.set_title(letters[5], loc="left", fontsize=11)
density_scatter(x, y, ax)
#ax.plot(x, y, '.', color=tableau20[1])
ax.plot([0, 40], [0, 40], '-', color='r')
ax.set_aspect('equal', adjustable='box')
ax.set_ylabel(name)
ax.set_xlabel("TREE-QMC-wh (n2)")
ax.set_xticks([0,  20, 40])
ax.set_yticks([0,  20, 40])

diff = y - x
n1 = len(numpy.where(diff == 0)[0])
n2 = len(numpy.where(diff > 0)[0])   # positive means TREE-QMC is better
n3 = len(numpy.where(diff < 0)[0])   # negative means other method is better
print("      Norm RF --> TIE: %d, TQMC better: %d, %s better: %d, TOTAL: %d" % (n1, n2, name, n3, n1+n2+n3))


# g) TREE-QMC-wh (n2) vs. ASTRAL/ASTER-wh - SPECIES TREE ERROR - RF
name = "ASTRAL-wh"
y = df.ASTERHxSERF.values * 100

ax = ax13
ax.set_title(letters[6], loc="left", fontsize=11)
density_scatter(x, y, ax)
#ax.plot(x, y, '.', color=tableau20[1])
ax.plot([0, 40], [0, 40], '-', color='r')
ax.set_aspect('equal', adjustable='box')
ax.set_ylabel(name)
ax.set_xlabel("TREE-QMC-wh (n2)")
ax.set_xticks([0,  20, 40])
ax.set_yticks([0,  20, 40])

diff = y - x
n1 = len(numpy.where(diff == 0)[0])
n2 = len(numpy.where(diff > 0)[0])   # positive means TREE-QMC is better
n3 = len(numpy.where(diff < 0)[0])   # negative means other method is better
print("      Norm RF --> TIE: %d, TQMC better: %d, %s better: %d, TOTAL: %d" % (n1, n2, name, n3, n1+n2+n3))


# QUARTET SCORE
x = df.TQMCwhn2xQS.values
y = df.ASTERHxQS.values

diff = y - x
n1 = len(numpy.where(diff == 0)[0])
n2 = len(numpy.where(diff < 0)[0])   # negative means TREE-QMC is better
n3 = len(numpy.where(diff > 0)[0])   # positive means other method is better
print("           QS --> TIE: %d, TQMC better: %d, %s better: %d, TOTAL: %d" % (n1, n2, name, n3, n1+n2+n3))

# h) TREE-QMC-wh (n2) vs. ASTRAL/ASTER-wh - LOCAL PP
x = df.TQMCwhn2xAVGLPP.values
y = df.ASTERHxAVGLPP.values

ax = ax14
ax.set_title(letters[7], loc="left", fontsize=11)
density_scatter(x, y, ax)
#ax.plot(x, y, '.', color=tableau20[1])
ax.plot([0.8, 1.0], [0.8, 1.0], '-', color='r')
ax.set_aspect('equal', adjustable='box')
ax.set_ylabel(name)
ax.set_xlabel("TREE-QMC-wh (n2)")
ax.set_xticks([0.8, 0.9, 1.0])
ax.set_yticks([0.8, 0.9, 1.0])

diff = y - x
n1 = len(numpy.where(diff == 0)[0])
n2 = len(numpy.where(diff < 0)[0])   # negative means TREE-QMC is better
n3 = len(numpy.where(diff > 0)[0])   # positive means other method is better
print("Mean Local PP --> TIE: %d, TQMC better: %d, %s better: %d, TOTAL: %d" % (n1, n2, name, n3, n1+n2+n3))

###########

# Mirarab et al (2015) 
# Read CSVs and remove sh support (focus on abayes), de-duplicate, and concatenate
print("\nMirarab et al (2015)")
df1 = pandas.read_csv("mirarab2015astral2-extsim/csvs/data-varyntax-for-testing.csv")
df2 = pandas.read_csv("mirarab2015astral2-extsim/csvs/data-varyils-for-testing.csv")
df2 = df2.drop(df2[(df2["NTAX"] == 200) &
                   (df2["ILSL"] == "medium") &
                   (df2["SPEC"] == "shallow")].index)
df = pandas.concat([df1, df2])

ax = ax20
ax.text(0.0, 0.5, "Data sets from\nMirarab et al. (2015)\nw/ abayes support", fontsize=12)
ax.spines["bottom"].set_visible(False)
ax.spines["left"].set_visible(False)
ax.set_xticks([])
ax.set_yticks([])

# i, j, k)
x = df.TQMCwhn2xSERF.values * 100

# i) TREE-QMC (wh, n2) vs. unweighted TREE-QMC (wf, n2)
name = "unweighted (n2)"
y = df.TQMCn2xSERF.values * 100

ax = ax21
ax.set_title(letters[8], loc="left", fontsize=11)
density_scatter(x, y, ax)
#ax.plot(x, y, '.', color=tableau20[1])
ax.plot([0, 40], [0, 40], '-', color='r')
ax.set_aspect('equal', adjustable='box')
ax.set_ylabel(name)
ax.set_xlabel("TREE-QMC-wh (n2)\n\n\% RF Error")
ax.set_xticks([0, 20, 40])
ax.set_yticks([0, 20, 40])

diff = y - x
n1 = len(numpy.where(diff == 0)[0])
n2 = len(numpy.where(diff > 0)[0])   # positive means TREE-QMC is better
n3 = len(numpy.where(diff < 0)[0])   # negative means other method is better
print("      Norm RF --> TIE: %d, TQMC better: %d, %s better: %d, TOTAL: %d" % (n1, n2, name, n3, n1+n2+n3))


# j) TREE-QMC (wh, n2) vs. ASTRID (ws)
name = "ASTRID-ws"
y = df.WASTRIDxSERF.values * 100

ax = ax22
ax.set_title(letters[9], loc="left", fontsize=11)
density_scatter(x, y, ax)
#ax.plot(x, y, '.', color=tableau20[1])
ax.plot([0, 40], [0, 40], '-', color='r')
ax.set_aspect('equal', adjustable='box')
ax.set_xlabel("TREE-QMC-wh (n2)\n\n\% RF Error")
ax.set_ylabel(name)
ax.set_xticks([0, 20, 40])
ax.set_yticks([0, 20, 40])

diff = y - x
n1 = len(numpy.where(diff == 0)[0])
n2 = len(numpy.where(diff > 0)[0])   # positive means TREE-QMC is better
n3 = len(numpy.where(diff < 0)[0])   # negative means other method is better
print("      Norm RF --> TIE: %d, TQMC better: %d, %s better: %d, TOTAL: %d" % (n1, n2, name, n3, n1+n2+n3))


name = "CA-ML"
y = df.CAMLxSERF.values * 100
diff = y - x
n1 = len(numpy.where(diff == 0)[0])
n2 = len(numpy.where(diff > 0)[0])   # positive means TREE-QMC is better
n3 = len(numpy.where(diff < 0)[0])   # negative means other method is better
print("      Norm RF --> TIE: %d, TQMC better: %d, %s better: %d, TOTAL: %d" % (n1, n2, name, n3, n1+n2+n3))


# k) TREE-QMC (wh, n2) vs. ASTER (wh)
name = "ASTER-wh"
y = df.ASTERHxSERF.values * 100

ax = ax23
ax.set_title(letters[10], loc="left", fontsize=11)
density_scatter(x, y, ax)
#ax.plot(x, y, '.', color=tableau20[1])
ax.plot([0, 40], [0, 40], '-', color='r')
ax.set_aspect('equal', adjustable='box')
ax.set_xlabel("TREE-QMC-wh (n2)\n\n\% RF Error")
ax.set_ylabel(name)
ax.set_xticks([0, 20, 40])
ax.set_yticks([0, 20, 40])

diff = y - x
n1 = len(numpy.where(diff == 0)[0])
n2 = len(numpy.where(diff > 0)[0])   # positive means TREE-QMC is better
n3 = len(numpy.where(diff < 0)[0])   # negative means other method is better
print("      Norm RF --> TIE: %d, TQMC better: %d, %s better: %d, TOTAL: %d" % (n1, n2, name, n3, n1+n2+n3))

# QUARTET SCORE
x = df.TQMCwhn2xQS.values
y = df.ASTERHxQS.values

diff = y - x
n1 = len(numpy.where(diff == 0)[0])
n2 = len(numpy.where(diff < 0)[0])   # negative means TREE-QMC is better
n3 = len(numpy.where(diff > 0)[0])   # positive means other method is better
print("           QS --> TIE: %d, TQMC better: %d, %s better: %d, TOTAL: %d" % (n1, n2, name, n3, n1+n2+n3))

# l) LOCAL PP
x = df.TQMCwhn2xAVGLPP.values
y = df.ASTERHxAVGLPP.values

ax = ax24
ax.set_title(letters[11], loc="left", fontsize=11)
density_scatter(x, y, ax)
#ax.plot(x, y, '.', color=tableau20[1])
ax.plot([0.8, 1.0], [0.8, 1.0], '-', color='r')
ax.set_aspect('equal', adjustable='box')
ax.set_xlabel("TREE-QMC-wh (n2)\n\nLocal PP")
ax.set_ylabel("ASTER-wh")
ax.set_xticks([0.8, 0.9, 1.0])
ax.set_yticks([0.8, 0.9, 1.0])

diff = y - x
n1 = len(numpy.where(diff == 0)[0])
n2 = len(numpy.where(diff < 0)[0])   # negative means TREE-QMC is better
n3 = len(numpy.where(diff > 0)[0])   # positive means other method is better
print("Mean Local PP --> TIE: %d, TQMC better: %d, %s better: %d, TOTAL: %d" % (n1, n2, name, n3, n1+n2+n3))

###########

for ax in axs:
    ax.tick_params(axis='x', labelsize=9)
    ax.tick_params(axis='y', labelsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


gs.tight_layout(fig, rect=[0, 0, 1, 1])
plt.savefig("plot-error-per-data-set.pdf", format='pdf', dpi=300)


"""
Morel et al (2023)
      Norm RF --> TIE: 84,  TQMC better: 1117,  shared better:   99, TOTAL: 1300
      Norm RF --> TIE: 110, TQMC better: 966,   ASTRID better:  224, TOTAL: 1300
      Norm RF --> TIE: 259, TQMC better: 442, ASTEROID better:  599, TOTAL: 1300
      Norm RF --> TIE: 178, TQMC better: 801,    ASTER better:  321, TOTAL: 1300
           QS --> TIE:  21, TQMC better:  13,    ASTER better: 1266, TOTAL: 1300
Mean Local PP --> TIE:   9, TQMC better: 725,    ASTER better:  566, TOTAL: 1300

# 4 nbps, 4 ngen, 50 replicates = 800 data sets
Zhang et al (2022)
      Norm RF --> TIE: 212, TQMC better: 453, unweighted better: 135, TOTAL: 800
      Norm RF --> TIE: 207, TQMC better: 421,     ASTRID better: 172, TOTAL: 800
      Norm RF --> TIE: 279, TQMC better: 351,      ASTER better: 170, TOTAL: 800
           QS --> TIE: 152, TQMC better:   0,      ASTER better: 648, TOTAL: 800
Mean Local PP --> TIE: 151, TQMC better:  41,      ASTER better: 608, TOTAL: 800

# (6 + 5 model conditions) * 3 ngen * 50 = 1650 - 31 (removed) = 1619 data sets
# CA-ML is minus an additional 48 for the 1000-taxon, 1000-gene condition (bc removed 2 for ASTER 1 thread)

Mirarab et al (2015)
      Norm RF --> TIE: 393, TQMC better: 941, unweighted better: 285, TOTAL: 1619
      Norm RF --> TIE: 474, TQMC better: 867,  ASTRID-ws better: 278, TOTAL: 1619
      Norm RF --> TIE: 269, TQMC better: 1038,     CA-ML better: 264, TOTAL: 1571
      Norm RF --> TIE: 597, TQMC better: 742,   ASTER-wh better: 280, TOTAL: 1619
           QS --> TIE: 392, TQMC better: 1,     ASTER-wh better: 1226, TOTAL: 1619
Mean Local PP --> TIE: 310, TQMC better: 73,    ASTER-wh better: 1236, TOTAL: 1619
"""

