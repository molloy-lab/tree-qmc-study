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
darkblue_pair = [(31, 119, 180), (174, 199, 232)]
orange_pair = [(255, 127, 14), (255, 187, 120)]
green_pair = [(44, 160, 44), (152, 223, 138)]
red_pair = [(214, 39, 40), (255, 152, 150)]
purple_pair = [(148, 103, 189), (197, 176, 213)]
brown_pair = [(140, 86, 75), (196, 156, 148)]
pink_pair = [(227, 119, 194), (247, 182, 210)]
gray_pair = [(127, 127, 127), (199, 199, 199)]
yellow_pair = [(188, 189, 34), (219, 219, 141)]
lightblue_pair = [(23, 190, 207), (158, 218, 229)]
specialgray = [(76.5, 76.5, 76.5)]

def map_rgb_to_01(tableau20):
    """
    Map RGB values to the [0, 1]
    """
    for i in range(len(tableau20)):
        r, g, b = tableau20[i]
        tableau20[i] = (r / 255., g / 255., b / 255.)

def setBoxColors(bp, tableau20, addfliers=False):
    """
    Color grouped boxplots
    """
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
        if addfliers:
            plt.setp(bp['fliers'][i], marker=".",
                     markerfacecolor=[0.3, 0.3, 0.3],
                     markeredgecolor=[0.3, 0.3, 0.3],
                     markersize=4)


def make_figure(df, supp, output):
    # Define colors for methods
    mthds = ["ASTRID-ws",
             "ASTER-wh",
             "TQMC-wh_n2",
             "TQMC-n2",
             "TQMC-wh_n0"]
    n_mthds = len(mthds)

    tableau20 = []
    tableau20 += orange_pair
    tableau20 += purple_pair
    tableau20 += darkblue_pair
    tableau20 += brown_pair
    tableau20 += pink_pair
    map_rgb_to_01(tableau20)

    # Create a figure
    fig = plt.figure(figsize=(7, 7.125))   # x,y
    gs = gridspec.GridSpec(3, 1)  # nrows, ncols
    ax00 = plt.subplot(gs[0, 0])
    ax10 = plt.subplot(gs[1, 0])
    ax20 = plt.subplot(gs[2, 0])
    axs = [ax00, ax10, ax20]

    ngens = [1000, 200, 50]        # Subplots are number of genes
    nbpss = [200, 400, 800, 1600]  # Box plots grouped by varying number of base pairs
    n_nbps = len(nbpss)

    for sub_ind, ngen in enumerate(ngens):
        ax = axs[sub_ind]

        # Set up data for box plot
        sers = [None] * n_nbps
        nrps = [None] * n_nbps
        for j, nbps in enumerate(nbpss):
            #print(nbps)
            sers[j] = []
            nrps[j] = []
            for k, mthd in enumerate(mthds):
                #print(mthd)
                if (mthd == "TQMC-n2"):
                    # IMPORTANT: TQMC-n2 doesn't use support
                    # and bs support means unrefined trees
                    # that are internally refined (better than abayes)!
                    ydf = df[(df["NBPS"] == nbps) &
                             (df["NGEN"] == ngen) &
                             (df["MTHD"] == mthd) &
                             (df["SUPP"] == "bs")]
                else:
                    ydf = df[(df["NBPS"] == nbps) &
                             (df["NGEN"] == ngen) &
                             (df["MTHD"] == mthd) &
                             (df["SUPP"] == supp)]

                ydf = ydf.sort_values(by=["REPL"], ascending=True)
                ser = ydf.SERF.values * 100
                sers[j].append(list(ser))
                nrps[j].append(len(ser))

        # Check data
        for j, nbps in enumerate(nbpss):
            if nrps[j][0] != 50:
                sys.stdout.write("%d genes, %d bp - Found %d replicates!\n" % (ngen, nbps, nrps[j][0]))

        # Create box plot
        xs = []
        ys = []
        inds = nbpss
        for ind, ser in zip(inds, sers):
            if ser != []:
                xs.append(ind)
                ys.append(ser)
        labs = xs
        sers = ys

        xminor = []
        xmajor = []
        base = numpy.arange(1, n_mthds + 1) 
        for j, nbps in enumerate(nbpss):
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

        if sub_ind == 0:
            mytitle = str("S100 Data (%s support)" % supp)
            #ax.set_title(mytitle, 
            #             loc="center", y=1.1, fontsize=12)
        elif sub_ind == 2:
            ax.set_xlabel("Sequence Length", fontsize=12)

        # Set tick labels
        ax.set_xlim(xminor[0]-1, xminor[-1]+1)
        ax.set_xticks(xmajor)
        test = []
        for j, nbps in enumerate(nbpss):
            test.append(r"%d" % nbps)
            #test.append(r"%d (%d)" % (nbps, nrps[j][0]))
        ax.set_xticklabels(test)

        # Set dashed lines
        if ngen == 1000:
            yticks = list(range(0, 20 + 1, 5))
        elif ngen == 200:
            yticks = list(range(0, 25 + 1, 5))
        else:
            yticks = list(range(0, 35 + 1, 5))

        # Add dashed lines
        ymax = yticks[-1]
        ymin = -0.05 * ymax
        ax.set_ylim(ymin, ymax)
        ax.set_yticks(yticks)
    
        xs = numpy.arange(xminor[0]-1, xminor[-1]+2)
        for y in yticks:
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
df = pandas.read_csv("../csvs/../csvs/data-all-error.csv")
title = str("plot-s100-varynbps-%s-supp.pdf" % (supp))
make_figure(df, supp, title)
