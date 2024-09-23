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

def make_figure(df, supp, output):
    mthds = ["ASTRID-ws",
             "ASTER-wh",
             "TQMC-wh_n2",
             "TQMC-n2",
             "CA-ML"]

    tableau20 = []
    tableau20 += orange
    tableau20 += purple
    tableau20 += darkblue
    tableau20 += brown
    tableau20 += red
    map_rgb_to_01(tableau20)

    n_mthds = len(mthds)

    fig = plt.figure(figsize=(8, 7.125))   # x,y
    gs = gridspec.GridSpec(3, 1)  # nrows, ncols
    ax00 = plt.subplot(gs[0, 0])  ## changing number of taxa
    ax10 = plt.subplot(gs[1, 0])
    ax20 = plt.subplot(gs[2, 0])
    axs = [ax00, ax10, ax20]

    # Plot error for varying number of taxa
    modls = ["high ILS\n(deep)",   "high ILS\n(shallow)",
             "medium ILS\n(deep)", "medium ILS\n(shallow)",
             "low ILS\n(deep)",    "low ILS\n(shallow)",]

    ilsls = ["high", "medium", "low"]
    specs = ["deep", "shallow"] 

    ngens = [1000, 200, 50]
    n_modl = len(modls)

    for sub_ind, ngen in enumerate(ngens):
        ax = axs[sub_ind]
        sers = [None] * n_modl
        nrps = [None] * n_modl
        j = 0
        for j1, ilsl in enumerate(ilsls):
            for j2, spec in enumerate(specs):
                #print(ntax)
                sers[j] = []
                nrps[j] = []
                for k, mthd in enumerate(mthds):
                    print(mthd)
                    ydf = df[(df["ILSL"] == ilsl) &
                             (df["SPEC"] == spec) &
                             (df["NGEN"] == ngen) &
                             (df["MTHD"] == mthd)]

                    ydf = ydf.sort_values(by=["REPL"], ascending=True)
                    #ser = ydf.SERF.values
                    ser = ydf.SEFNR.values * 100
                    sers[j].append(list(ser))
                    nrps[j].append(len(ser))
                j += 1

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
        base = numpy.arange(1, n_mthds + 1) 
        for j, modl in enumerate(modls):
            pos = base + ((n_mthds + 1) * j)
            xminor = xminor + list(pos)
            xmajor = xmajor + [numpy.mean(pos)]

            bp = ax.boxplot(sers[j], positions=pos, widths=0.75,
                            showfliers=False, 
                            showmeans=True,
                            patch_artist=True)
            setBoxColors(bp, tableau20)

        # Set labels
        ax.set_title(letters[sub_ind] + str("   %d genes" % ngen),
                     loc="left", x=0.0, y=1.0, fontsize=10.5)
        ax.set_ylabel(r"\% RF Error", fontsize=11)  

        ax.tick_params(axis='x', labelsize=10)
        ax.tick_params(axis='y', labelsize=9)        

        # Set tick labels
        ax.set_xlim(xminor[0]-1, xminor[-1]+1)
        ax.set_xticks(xmajor)
        test = []
        for j, modl in enumerate(modls):
            if nrps[j][0] != 50:
                sys.stdout.write("modl %s, ngen %d - Found %d replicates!\n" % (modl, ngen, nrps[j][0]))
            test.append(r"%s" % modl)

        ax.set_xticklabels(test)

        # Set dashed lines
        if ngen == 1000:
            ymin = 0
            ymax = 15
            yticks = list(range(0, ymax + 1, 2))
        elif ngen == 200:
            ymin = 0
            ymax = 30
            yticks = list(range(0, ymax + 1, 5))
        else:
            ymin = 0
            ymax = 40
            yticks = list(range(0, ymax + 1, 5))
        
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

    # Shift layout
    gs.tight_layout(fig, rect=[0, 0.05, 1, 1])

    # Add legend at bottom
    hs = []
    for k in range(n_mthds):
        #print(mthds[k])
        h, = ax.plot([1], [1], 
                     '-', 
                     color=tableau20[k*2], 
                     lw=10)

        hs.append(h)
        ax.legend(hs, mthds,
                  frameon=False, 
                  ncol=7, 
                  fontsize=9.5,
                  loc='lower center', 
                  bbox_to_anchor=(0.5, -0.5, 0, 3))

    # Save plot
    plt.savefig(output, format='pdf', dpi=300)

# Read and plot data
supp = "abayes"
df = pandas.read_csv("../../csvs/data-varyils-error-and-qscore.csv")
title = str("plot-astral2-varyils-%s-supp.pdf" % (supp))
make_figure(df, supp, title)




