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

fig = plt.figure(figsize=(5, 2))  # x, y
gs = gridspec.GridSpec(1, 2)       # nrows, ncols (reversed)

ax00 = plt.subplot(gs[0, 0])
ax01 = plt.subplot(gs[0, 1])

axs = [ax00, ax01]

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
x = df.TQMCn2refinedxFNR.values * 100  # refined i.e. binary so FN rate = RF rate

# Only plot TREE-QMC (n2) vs. Asteroid
name = "Asteroid"
y = df.ASTEROIDxFNR.values * 100

ax = ax01
ax.set_title(letters[2], loc="left", fontsize=11)
density_scatter(x, y, ax)
#ax.plot(x, y, '.', color=tableau20[1])
ax.plot([0, 100], [0, 100], '-', color='r')
ax.set_aspect('equal', adjustable='box')
ax.set_ylabel(name)
ax.set_xlabel("TREE-QMC (n2)")
ax.set_xticks([0,  50, 100])
ax.set_yticks([0,  50, 100])

###########

for ax in axs:
    ax.tick_params(axis='x', labelsize=9)
    ax.tick_params(axis='y', labelsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


gs.tight_layout(fig, rect=[0, 0, 1, 1])
plt.savefig("plot-error-morel-data-set-asteroid.pdf", format='pdf', dpi=300)


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

