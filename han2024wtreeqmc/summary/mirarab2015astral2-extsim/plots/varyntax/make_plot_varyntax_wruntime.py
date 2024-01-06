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
    ymin = yticks[0]
    ymax = yticks[-1]

    ymin = -0.05 * ymax
    ax.set_ylim(ymin, ymax)
    ax.set_yticks(yticks)

    if xminor == None:
        xdraw = ax.xaxis.get_majorticklocs()
        for y in yticks:
            ax.plot(xdraw, [y] * len(xdraw), "--", dashes=(3, 3),
                    lw=0.5, color="black", alpha=0.3)
    else:
        xs = numpy.arange(xminor[0]-1, xminor[-1]+2)
        for y in yticks:
            ax.plot(xs, [y] * len(xs), "--", dashes=(3, 3),
                    lw=0.5, color="black", alpha=0.3)


def plot_error(ax, df, ntaxs, supp):
    n_ntax = len(ntaxs)

    ax.set_title(letters[0],
                 loc="left", x=0.0, y=1.0, fontsize=10.5)
    #ax.set_title(upperletters[0],
    #             loc="left", fontsize=11)

    mthds = ["ASTRID-ws",
             "ASTER-wh",
             "TQMC-wh_n2",
             "TQMC-n2",
             "CA-ML"]
    n_mthds = len(mthds)

    tableau20 = []
    tableau20 += orange
    tableau20 += purple
    tableau20 += darkblue
    tableau20 += brown
    tableau20 += red

    map_rgb_to_01(tableau20)

    sers = [None] * n_ntax
    nrps = [None] * n_ntax
    for j, ntax in enumerate(ntaxs):
        #print(ntax)
        sers[j] = []
        nrps[j] = []
        for k, mthd in enumerate(mthds):
            print(mthd)
            if (mthd == "CA-ML") or (mthd == "TQMC-n2"):
                ydf = df[(df["NTAX"] == ntax) &
                          (df["NGEN"] == 1000) &
                          (df["MTHD"] == mthd) &
                          (df["SUPP"] == "none")]
            else: 
                ydf = df[(df["NTAX"] == ntax) &
                          (df["NGEN"] == 1000) &
                          (df["MTHD"] == mthd) &
                          (df["SUPP"] == supp)]

            ydf = ydf.sort_values(by=["REPL"], ascending=True)
            #ser = ydf.SERF.values
            ser = ydf.SEFNR.values * 100
            #print(ser)
            sers[j].append(list(ser))
            nrps[j].append(len(ser))

    xs = []
    ys = []
    inds = ntaxs
    for ind, ser in zip(inds, sers):
        if ser != []:
            xs.append(ind)
            ys.append(ser)
    labs = xs
    sers = ys

    xminor = []
    xmajor = []
    base = numpy.arange(1, n_mthds + 1) 
    for j, ntax in enumerate(ntaxs):
        pos = base + ((n_mthds + 1) * j)
        xminor = xminor + list(pos)
        xmajor = xmajor + [numpy.mean(pos)]

        bp = ax.boxplot(sers[j], positions=pos, widths=0.75,
                        showfliers=False, 
                        showmeans=True,
                        patch_artist=True)
        setBoxColors(bp, tableau20)

    # Set labels
    ax.set_ylabel(r"\% RF Error", fontsize=11)  
    ax.tick_params(axis='x', labelsize=10)
    ax.tick_params(axis='y', labelsize=9)
    
    # Set tick labels
    ax.set_xlim(xminor[0]-1, xminor[-1]+1)
    ax.set_xticks(xmajor)
    test = []
    for j, ntax in enumerate(ntaxs):
        if nrps[j][0] != 50:
            sys.stdout.write("ntax %d - Found %d replicates!\n" % (ntax, nrps[j][0]))
        test.append(r"%d" % ntax)

    ax.set_xticklabels(test)

    # Set dashed lines
    yticks = list(range(0, 10 + 1, 2))
    add_dashed_lines(ax, yticks, xminor=xminor)

    # Set plot axis parameters
    ax.tick_params(axis=u'both', which=u'both',length=0) # removes tiny ticks
    ax.get_xaxis().tick_bottom() 
    ax.get_yaxis().tick_left()

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    return [mthds, tableau20]


def plot_runtime(ax, df, ntaxs, supp):
    # Plot runtime for number of taxa
    n_ntax = len(ntaxs)

    ax.set_title(letters[1],
                 loc="left", x=0.0, y=1.0, fontsize=10.5)
    #ax.set_title(upperletters[1],
    #             loc="left", fontsize=11)

    mthds = ["ASTRID-ws",
             "ASTER-wh (1 thread)",
             "ASTER-wh (16 threads)",
             "TQMC-wh_n2"]
    n_mthds = len(mthds)

    tableau20 = []
    tableau20 += orange
    tableau20 += purple
    tableau20 += purple
    tableau20 += darkblue
    map_rgb_to_01(tableau20)

    for k, mthd in enumerate(mthds):
        print(mthd)
        av = []
        se = []
        for ntax in ntaxs:
            xdf = df2[(df2["NTAX"] == ntax) & 
                      (df2["NGEN"] == 1000) &
                      (df2["SUPP"] == "sh") &  # update after jobs are done
                      (df2["MTHD"] == mthd)]

            if xdf.shape[0] > 0:
                vals = xdf.SECS.values / (60.0 * 60.0)

                av.append(numpy.mean(vals))
                se.append(numpy.std(vals) / numpy.sqrt(vals.size))
    
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

    yticks = range(0, 9, 2)
    add_dashed_lines(ax, yticks)

    ax.set_xlabel("Number of Taxa", fontsize=12)
    ax.set_ylabel(r'Runtime (h)', fontsize=11)
    ax.tick_params(axis='x', labelsize=10)
    ax.tick_params(axis='y', labelsize=9)

    ax.get_xaxis().tick_bottom() 
    ax.get_yaxis().tick_left() 
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    #print("legend")
    #h1, = ax.plot([1, 2, 3], [1, 2, 3], '-', color='k', lw=1)
    #h2, = ax.plot([1, 2, 3], [1, 2, 3], '--', color='k', lw=1)
    #ax.legend([h1, h2],
    #          ["1 thread", "16 threads"], 
    #          frameon=False,
    #          ncol=2,
    #          fontsize=8,
    #          loc='upper left')

    ax.text(0, 7, 
            "Solid lines = 1 thread\nDashed line = 16 threads",
             fontsize=9,
             horizontalalignment='left',
             verticalalignment='center')


def make_figure(df1, df2, supp, output):
    fig = plt.figure(figsize=(8, 4.75))   # x,y
    gs = gridspec.GridSpec(2, 1)  # nrows, ncols
    ax00 = plt.subplot(gs[0, 0])  ## changing number of taxa
    ax10 = plt.subplot(gs[1, 0])

    ntaxs = [10, 50, 100, 200, 500, 1000]

    plot_runtime(ax10, df2, ntaxs, supp)
    [mthds, tableau20] = plot_error(ax00, df1, ntaxs, supp)

    # Shift layout
    gs.tight_layout(fig, rect=[0, 0.05, 1, 1])

    # Add legend at bottom
    hs = []
    for k in range(len(mthds)):
        #print(mthds[k])
        h, = ax00.plot([0], [0], 
                       '-', 
                       color=tableau20[k*2], 
                       lw=10)
        hs.append(h)

    ax10.legend(hs, mthds,
                frameon=False,
                ncol=5, 
                fontsize=9.5,
                loc='lower center', 
                bbox_to_anchor=(0.5, -0.575, 0, 3))

    # Save plot
    plt.savefig(output, format='pdf', dpi=300)


# Read and plot data
supp = "abayes"
df1 = pandas.read_csv("../../csvs/data-varyntax-error.csv",
                      keep_default_na=True)
df2 = pandas.read_csv("../../csvs/data-varyntax-error-and-runtime.csv",
                      keep_default_na=True, quoting=csv.QUOTE_NONNUMERIC)
title = str("plot-astral2-varyntax-%s-supp-wruntime.pdf" % (supp))
make_figure(df1, df2, supp, title)
