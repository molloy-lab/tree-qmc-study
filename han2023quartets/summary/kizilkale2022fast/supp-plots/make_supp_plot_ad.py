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
           r"\textbf{f)}",
           r"\textbf{g)}", 
           r"\textbf{h)}", 
           r"\textbf{i)}"]

upperletters = [r"\textbf{A}", 
                r"\textbf{B}", 
                r"\textbf{C}", 
                r"\textbf{D}", 
                r"\textbf{E}", 
                r"\textbf{F}",
                r"\textbf{G}", 
                r"\textbf{H}", 
                r"\textbf{I}",]

# Tableau 20 colors in RGB.    
tableau20 = [(44, 160, 44), (152, 223, 138),
             (255, 127, 14), (255, 187, 120),
             (23, 190, 207), (158, 218, 229),
             (31, 119, 180), (174, 199, 232),
             (227, 119, 194), (247, 182, 210),
             (214, 39, 40), (255, 152, 150),
             (148, 103, 189), (197, 176, 213)]

# Dark blue - (31, 119, 180), (174, 199, 232),
# Orange - (255, 127, 14), (255, 187, 120),
# Green - (44, 160, 44), (152, 223, 138),
# Red - (214, 39, 40), (255, 152, 150), 
# Purple - (148, 103, 189), (197, 176, 213),
# Brown - (140, 86, 75), (196, 156, 148), 
# Pink - (227, 119, 194), (247, 182, 210),
# Gray - (127, 127, 127), (199, 199, 199), 
# Yellow - (188, 189, 34), (219, 219, 141),
# Light blue - (23, 190, 207), (158, 218, 229)

# Map RGB values to the [0, 1]
for i in range(len(tableau20)):
    r, g, b = tableau20[i]
    tableau20[i] = (r / 255., g / 255., b / 255.)


# Color grouped boxplots
def setBoxColors(bp):
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
        plt.setp(bp['fliers'][i], marker=".",
                 markerfacecolor=[0.3, 0.3, 0.3],
                 markeredgecolor=[0.3, 0.3, 0.3],
                 markersize=4)


# https://stackoverflow.com/questions/25812255/row-and-column-headers-in-matplotlibs-subplots
# https://stackoverflow.com/questions/27426668/row-titles-for-matplotlib-subplot


def make_figure(df, output):
    mthds = ["huntress_v0.1.2.0_default",
             "fastral_wrootx_wmuts",
             "treeqmcbip_v1.0.0_n2_wrootx_wmuts",
             "scistree_v1.2.0.6_wmuts",
             "fastme_v2.1.5_wrootx_wmuts"]
    mnams = ["HUNTRESS",
             r"FASTRAL",
             r"TREE-QMC",
             r"ScisTree",
             r"FastME"]
    m = len(mthds)

    ncxnms = [["n1000", "m300"], ["n300", "m300"], ["n300", "m1000"]]
    n = len(ncxnms)

    #fig = plt.figure(figsize=(8, 4.75))
    fig = plt.figure(figsize=(6, 4))

    gs = gridspec.GridSpec(3,3)

    ax00 = plt.subplot(gs[0,0])  ## precision  300 x 300
    ax01 = plt.subplot(gs[0,1])  ## precision  300 x 1000
    ax02 = plt.subplot(gs[0,2])  ## precision 1000 x 300

    ax10 = plt.subplot(gs[1,0])  ## recall  300 x 300
    ax11 = plt.subplot(gs[1,1])  ## recall  300 x 1000
    ax12 = plt.subplot(gs[1,2])  ## recall 1000 x 300

    ax20 = plt.subplot(gs[2,0])  ## flip  300 x 300
    ax21 = plt.subplot(gs[2,1])  ## flip  300 x 1000
    ax22 = plt.subplot(gs[2,2])  ## flip 1000 x 300

    grid = [[ax00, ax01, ax02],
            [ax10, ax11, ax12],
            [ax20, ax21, ax22]]

    for i in range(3):
        for j in range(n):
            print("Ax %d %d" % (i, j))
            ax = grid[i][j]

            ncell = ncxnms[j][0]
            nmuts = ncxnms[j][1]
            print("  %s x %s" % (ncell, nmuts))
        
            sers = []
            avgs = []
            stds = []
            nrps = []
            for k, mthd in enumerate(mthds):
                ydf = df[(df["NCELL"] == ncell) &
                         (df["NMUT"] == nmuts) &
                         (df["BETA_FN"] == "fn0.2") &
                         (df["ALPHA_FP"] == "fp0.001") &
                         (df["GAMMA_NA"] == "na0.05") &
                         (df["S"] == "s100") &
                         (df["H"] == "h1") &
                         (df["MINVAF"] == "minVAF0.005") &
                         (df["ISAV"] == "ISAV0") &
                         (df["D"] == "d0") &
                         (df["L"] == "l1000000") &
                         (df["MTHD"] == mthd)]

                if i == 0:
                    vals = ydf["AD_P"].values
                elif i == 1:
                    vals = ydf["AD_R"].values
                else:
                    vals = ydf["AD_FLIP"].values
 
                sers.append(list(vals))
                avgs.append(numpy.mean(vals))
                stds.append(numpy.mean(vals))
                nrps.append(len(vals))

            pos = numpy.arange(1, m+1)
            if i == 2:
                # Bar graph
                poss = [2, 4, 6, 8, 10]
                for k, mthd in enumerate(mthds):
                    ax.bar([poss[k]],
                           [avgs[k]],
                           1.5,
                           yerr=[[0],[stds[k]]],
                           color=[tableau20[2*k+1]],
                           edgecolor=[tableau20[2*k+0]],
                           lw=1.5,
                           error_kw=dict(ecolor=tableau20[2*k+0],
                           capsize=3,
                           capthick=1.5))
            else:
                # Box plot
                bp = ax.boxplot(sers, positions=pos, widths=0.75,
                                showfliers=True, 
                                showmeans=True,
                                patch_artist=True)
                setBoxColors(bp)

            # Add title
            if i == 0:
                [ncell, nmuts] = ncxnms[j]
                ncell = ncell.replace('n', '')
                nmuts = nmuts.replace('m', '')
                labl = ncell + " cells x " + nmuts + " muts"
                ax.set_title(labl,
                             loc="center", x=0.5, y=1.3,
                             fontsize=11)

            # Add letter
            ax.text(0.05, 1.125, letters[i*3 + j], fontsize=10, 
                    horizontalalignment='center',
                    verticalalignment='center',
                    transform=ax.transAxes)

            # Set x ticks
            ax.set_xticklabels(['', '', '', '', ''], fontsize=9)

            if j == 0:
                if i == 0:
                    ax.set_ylabel("AD Precision", fontsize=10)
                elif i == 1:
                    ax.set_ylabel("AD Recall", fontsize=10)
                else:
                    ax.set_ylabel("AD flip rate", fontsize=10)

            # Set y ticks
            if i == 0:
                ytick_min = 0.8
                ytick_max = 1.00
                yticks = [0.8, 0.85, 0.9, 0.95, 1.0]
            elif i == 1:
                ytick_min = 0.5
                ytick_max = 1.0
                yticks = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
            else:
                if j == 0:
                    ytick_min = 0.0
                    ytick_max = 0.01
                    yticks = [0.0, 0.005, 0.01]
                elif j == 1:
                    ytick_min = 0.0
                    ytick_max = 0.025
                    yticks = [0.0, 0.0125, 0.025]
                else:
                    ytick_min = 0.0
                    ytick_max = 0.008
                    yticks = [0.0, 0.002, 0.004, 0.006, 0.008]

            diff = ytick_max - ytick_min
            ymin = ytick_min - diff * 0.05
            ymax = ytick_max + diff * 0.05
            ax.set_ylim(ymin, ymax)
            ax.set_yticks(yticks)
            ax.tick_params(axis='y', labelsize=9)
    
            #xs = numpy.arange(xminor[0]-1, xminor[-1]+2)
            #for y in ydraw:
            #    ax.plot(xs, [y] * len(xs), "--", dashes=(3, 3),
            #            lw=0.5, color="black", alpha=0.3)

            #ax.tick_params(axis='y', labelsize=8)

            # Set plot axis parameters
            ax.tick_params(axis=u'both', which=u'both',length=0) # removes tiny ticks
            ax.get_xaxis().tick_bottom() 
            ax.get_yaxis().tick_left()

            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)

    # Add legend at bottom
    gs.tight_layout(fig, rect=[0, 0.07, 1, 1])
    hs = []
    for k in range(len(mnams)):
        print(mnams[k])
        h, = ax.plot([1], [1], 
                     '-', 
                     color=tableau20[k*2], 
                     lw=10)
        hs.append(h)

    ax.legend(hs, 
              mnams,
              frameon=False, 
              ncol=5, 
              fontsize=9,
              loc='lower center', 
              bbox_to_anchor=(-0.95, -0.6, 0, 1))

    # Save plot
    plt.savefig(output, format='pdf', dpi=300)

# Read and plot data
df = pandas.read_csv("../csvs/data-all-error-and-timings_exfig4-6.csv")

title = "supp-plot-ancestor-descendant.pdf"
make_figure(df, title)
