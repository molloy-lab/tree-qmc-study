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

def make_figure(df, output):
    mthds = ["TREE-QMC-wf_n0",
             "TREE-QMC-wf_n1_shared",
             "TREE-QMC-wf_n2",
             "Asteroid",
             "ASTER",
             "ASTRID"]

    names = ["TQMC-n0",
             "TQMC-n1 (shared)",
             "TQMC-n2",
             "Asteroid",
             "ASTER",
             "ASTRID"]

    tableau20 = []
    tableau20 += brown
    tableau20 += lightblue
    tableau20 += darkblue
    tableau20 += pink
    tableau20 += purple
    tableau20 += orange

    map_rgb_to_01(tableau20)

    m = len(mthds)

    modls = [0.05, 0.1, 1, 10, 100, 200] # varying branch length scaler

    fig = plt.figure(figsize=(8, 4.75))  # x,y
    gs = gridspec.GridSpec(2, 1)  # nrows, ncols
    ax00 = plt.subplot(gs[0, 0])  # changing ntax FN
    ax10 = plt.subplot(gs[1, 0])  # changing ntax FP

    for doit in [0, 1]:
        if doit == 0:
            ax = ax00
        else:
            ax = ax10

        n = len(modls)

        sers = [None] * n
        nrps = [None] * n
        for j, modl in enumerate(modls):
            sers[j] = []
            nrps[j] = []
            for k, mthd in enumerate(mthds):
                ydf = df[(df["BLSC"] == modl) &
                         (df["MTHD"] == mthd)]
                ydf = ydf.sort_values(by=["REPL"], ascending=True)

                if doit == 0:
                    ser = ydf.SEFNR.values
                else:
                    ser = ydf.SEFPR.values
                sers[j].append(list(ser))
                nrps[j].append(len(ser))

        xs = []
        ys = []
        inds = modls
        for ind, ser in zip(inds, sers):
            if ser != []:
                xs.append(ind)
                ys.append(ser)
        labs = xs
        sers = ys

        xminor = []
        xmajor = []
        base = numpy.arange(1, m + 1) 
        for j, modl in enumerate(modls):
            pos = base + ((m + 1) * j)
            xminor = xminor + list(pos)
            xmajor = xmajor + [numpy.mean(pos)]

            bp = ax.boxplot(sers[j],
                            positions=pos,
                            widths=0.75,
                            showfliers=False,
                            showmeans=True,
                            patch_artist=True)
            setBoxColors(bp, tableau20)

        # Set labels
        ax.set_title(letters[doit],
                 loc="left", fontsize=11)
        if doit == 0:
            ax.set_ylabel(r"FN Error", fontsize=11)  
        else:
            ax.set_ylabel(r"FP Error", fontsize=11)
            ax.set_xlabel(r"Branch Length Scaler", fontsize=12)

        ax.tick_params(axis='x', labelsize=10)
        ax.tick_params(axis='y', labelsize=9)
        ax.set_xlim(xminor[0]-1, xminor[-1]+1)
        ax.set_xticks(xmajor)

        test = []
        for j, modl in enumerate(modls):
            if nrps[j][0] != 50:
                sys.stdout.write("%s - Found %d replicates!\n" % (modl, nrps[j][0]))
            test.append(r"%s" % modl)
        ax.set_xticklabels(test)

        if doit == 0:
            yticks = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
        else:
            yticks = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
        ydraw = yticks
        ymax = yticks[-1]

        ymin = -0.05 * ymax
        ax.set_ylim(ymin, ymax)
        ax.set_yticks(yticks)
    
        xs = numpy.arange(xminor[0]-1, xminor[-1]+2)
        for y in ydraw:
            ax.plot(xs, [y] * len(xs), "--", dashes=(3, 3),
                    lw=0.5, color="black", alpha=0.3)

        # Set plot axis parameters
        ax.tick_params(axis=u'both', which=u'both',length=0) # removes tiny ticks
        ax.get_xaxis().tick_bottom() 
        ax.get_yaxis().tick_left()

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    if doit == 1:
        # Add legend at bottom
        gs.tight_layout(fig, rect=[0, 0.05, 1, 1])
        hs = []
        for k in range(len(names)):
            print(names[k])
            h, = ax.plot([1], [1], 
                         '-', 
                         color=tableau20[k*2], 
                         lw=10)
        
            hs.append(h)
            ax.legend(hs, names,
                      frameon=False,
                      ncol=7, 
                      fontsize=9,
                      loc='upper center', 
                      bbox_to_anchor=(0.5, -0.275))

    # Save plot
    #gs.tight_layout(fig, rect=[0, 0, 1, 1])
    plt.savefig(output, format='pdf', dpi=300)

# Read and plot data
df = pandas.read_csv("../../csvs/data-varyblsc-error-and-timings.csv")
make_figure(df, "plot-asteroid-varyblsc.pdf")
