import csv
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

upperletters = [r"\textbf{A}", 
                r"\textbf{B}", 
                r"\textbf{C}", 
                r"\textbf{D}", 
                r"\textbf{E}", 
                r"\textbf{F}"]

# Tableau 20 colors in RGB.    
#tableau20 = [(44, 160, 44), (152, 223, 138),
#             (255, 127, 14), (255, 187, 120),
#             (23, 190, 207), (158, 218, 229),
#             (31, 119, 180), (174, 199, 232),
#             (227, 119, 194), (247, 182, 210),
#             (214, 39, 40), (255, 152, 150),
#             (148, 103, 189), (197, 176, 213)]

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

# Color grouped boxplots
def setBoxColors(bp, tableau20):
    n = len(bp['boxes'])
    for i in range(n):
        j = i

        plt.setp(bp['boxes'][i], 
                 color=tableau20[j*2],
                 linewidth=1.75)
        plt.setp(bp['boxes'][i], 
                 facecolor=tableau20[j*2+1])
        plt.setp(bp['caps'][i*2:i*2+2],
                 color=tableau20[j*2],
                 linewidth=1.75)
        plt.setp(bp['whiskers'][i*2:i*2+2],
                 color=tableau20[j*2],
                 linewidth=1.75,
                 linestyle='-')
        plt.setp(bp['medians'][i],
                 color=[0.3, 0.3, 0.3],
                 linewidth=1.75)
                 #color=tableau20[i*2], linewidth=1.75)
        plt.setp(bp['means'][i],
                 markerfacecolor=[0.3, 0.3, 0.3],
                 markeredgecolor=[0.3, 0.3, 0.3],
                 markersize=3)
        #plt.setp(bp['fliers'][i], marker=".",
        #         markerfacecolor=[0.3, 0.3, 0.3],
        #         markeredgecolor=[0.3, 0.3, 0.3],
        #         markersize=4)


def add_dashed_lines(ax, yticks, xminor=None):
    if xminor == None:
        xdraw = ax.xaxis.get_majorticklocs()
        for y in yticks:
            ax.plot(xdraw, [y] * len(xdraw), "--", dashes=(3, 3),
                    lw=1, color="black", alpha=0.3)
    else:
        xs = numpy.arange(xminor[0]-1, xminor[-1]+2)
        for y in yticks:
            ax.plot(xs, [y] * len(xs), "--", dashes=(3, 3),
                    lw=1, color="black", alpha=0.3)


def plot_runtime(ax, df, ntaxs):
    # Plot runtime for number of taxa
    n_ntax = len(ntaxs)

    ax.set_title(letters[1],
                 loc="left", x=0.0, y=1.0, fontsize=10.5)
    #ax.set_title(upperletters[1],
    #             loc="left", fontsize=11)

    mthds = ["ASTER-wh (1 thread)",
             "ASTER-wh (16 threads)",
             "TQMC-wh_n2",
             "TQMC-n2",
             "ASTRID-ws"]
    n_mthds = len(mthds)

    tableau20 = []
    tableau20 += purple
    tableau20 += purple
    tableau20 += darkblue
    tableau20 += brown
    tableau20 += orange
    map_rgb_to_01(tableau20)

    for k, mthd in enumerate(mthds):
        print(mthd)
        av = []
        se = []
        for ntax in ntaxs:
            xdf = df[(df["NTAX"] == ntax) & 
                     (df["MTHD"] == mthd)]

            vals = xdf.SECS.values / (60.0 * 60.0)

            av.append(numpy.mean(vals))
            se.append(numpy.std(vals)) # / numpy.sqrt(vals.size))
    
        print(av)

        av = numpy.array(av)
        se = numpy.array(se)

        labls = ntaxs
        if mthd == "ASTER-wh (16 threads)":
            ax.plot(labls, av, '--', color=tableau20[2*k], lw=3, 
                    label=mthds[k])
        else:
            ax.plot(labls, av, '-', color=tableau20[2*k], lw=2, 
                    label=mthds[k])
        
        ax.fill_between(labls, av - se, av + se, 
                        color=tableau20[2*k+1], alpha=0.5)

    ax.set_xticks(ntaxs)
    ax.set_xticklabels([str(x) for x in ntaxs])

    yticks = [1, 2, 3, 4, 5,  8]
    ax.set_ylim(0, yticks[-1])
    ax.set_yticks([0] + yticks)
    ax.set_yticklabels(['0.0', '1.0', '2.0', '3.0', '4.0', '5.0',  '8.0'])
    add_dashed_lines(ax, yticks)

    ax.set_xlabel("Number of Taxa", fontsize=12)
    ax.set_ylabel(r'Runtime (h)', fontsize=12)
    ax.tick_params(axis='x', labelsize=10.5)
    ax.tick_params(axis='y', labelsize=10.5)

    ax.get_xaxis().tick_bottom() 
    ax.get_yaxis().tick_left() 
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.text(0, 6.5, 
            "Solid lines = 1 thread\nDashed line = 16 threads",
            fontsize=12,
            horizontalalignment='left',
            verticalalignment='center')

def plot_runtime_ratio(ax, df, ntaxs):
    # Plot runtime for number of taxa
    n_ntax = len(ntaxs)

    ax.set_title(letters[0],
                 loc="left", x=0.0, y=1.0, fontsize=10.5)
    #ax.set_title(upperletters[1],
    #             loc="left", fontsize=11)

    tableau20 = []
    tableau20 += darkblue
    map_rgb_to_01(tableau20)

    av = []
    se = []
    for ntax in ntaxs:
        xdf = df[(df["NTAX"] == ntax) & (df["MTHD"] == "TQMC-n2")]
        xdf = xdf.sort_values(by=["REPL"], ascending=True)
        xvals = xdf.SECS.values

        ydf = df[(df["NTAX"] == ntax) & (df["MTHD"] == "TQMC-wh_n2")]
        ydf = ydf.sort_values(by=["REPL"], ascending=True)
        yvals = ydf.SECS.values

        if ((xdf.REPL.values == ydf.REPL.values).all()):
            pass
        else:
            sys.exit("Cannot do division!")

        vals = yvals / xvals

        av.append(numpy.mean(vals))
        se.append(numpy.std(vals)) # / numpy.sqrt(vals.size))

    av = numpy.array(av)
    se = numpy.array(se)

    labls = ntaxs
    
    ax.fill_between(labls, av - se, av + se, 
                    color=tableau20[1], alpha=0.5)

    ax.scatter(labls, av, color=tableau20[0]) #, size=3)

    ax.set_xticks(ntaxs)
    ax.set_xticklabels([str(x) for x in ntaxs])

    yticks = [1.5, 2, 2.5, 3, 3.5, 4]
    ax.set_ylim(1, yticks[-1])
    ax.set_yticks([1] + yticks)
    add_dashed_lines(ax, yticks)

    #ax.set_xlabel("Number of Taxa", fontsize=12)
    ax.set_ylabel(r'TQMC runtime ratio', fontsize=12)
    ax.tick_params(axis='x', labelsize=10.5)
    ax.tick_params(axis='y', labelsize=10.5)

    ax.get_xaxis().tick_bottom() 
    ax.get_yaxis().tick_left() 
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def make_figure(df, output):
    fig = plt.figure(figsize=(8, 4.75))   # x,y
    gs = gridspec.GridSpec(2, 1)  # nrows, ncols
    ax00 = plt.subplot(gs[0, 0])  # changing number of taxa
    ax10 = plt.subplot(gs[1, 0])

    ntaxs = [10, 50, 100, 200, 500, 1000]

    plot_runtime_ratio(ax00, df, ntaxs)
    plot_runtime(ax10, df, ntaxs)

    # Shift layout
    gs.tight_layout(fig, rect=[0, 0.05, 1, 1])

    # Add legend at bottom
    mthds = ["ASTER-wh",
             "TREE-QMC-wh (n2)",
             "TREE-QMC (n2)",
             "ASTRID-ws"]

    tableau20 = []
    tableau20 += purple
    tableau20 += darkblue
    tableau20 += brown
    tableau20 += orange
    map_rgb_to_01(tableau20)

    hs = []
    for k in range(len(mthds)):
        #print(mthds[k])
        h, = ax10.plot([0], [0], 
                       '-', 
                       color=tableau20[k*2], 
                       lw=10)
        hs.append(h)

    ax10.legend(hs, mthds,
                frameon=False,
                ncol=5, 
                fontsize=11,
                loc='lower center', 
                bbox_to_anchor=(0.5, -0.55, 0, 3))

    # Save plot
    plt.savefig(output, format='pdf', dpi=300)


# Read and plot data
supp = "abayes"
df = pandas.read_csv("../../csvs/data-varyntax-runtime.csv", keep_default_na=True)
df = df.drop(df[(df["NGEN"] == 200)].index)
df = df.drop(df[(df["NGEN"] == 50)].index)
df = df.drop(df[(df["SUPP"] == "sh")].index)
title = str("plot-astral2-varyntax-%s-runtime.pdf" % (supp))
make_figure(df, title)

