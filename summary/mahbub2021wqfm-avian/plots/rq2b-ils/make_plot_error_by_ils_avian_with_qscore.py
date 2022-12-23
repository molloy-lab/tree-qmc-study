import numpy
import pandas
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import rc
import sys


plt.rc('text', usetex=True)
plt.rcParams['text.latex.preamble'] = [r'\usepackage{helvet} \usepackage{sfmath}']

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
tableau20 = [(255, 127, 14), (255, 187, 120),
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
        #plt.setp(bp['fliers'][i], marker=".",
        #         markerfacecolor=[0.3, 0.3, 0.3],
        #         markeredgecolor=[0.3, 0.3, 0.3],
        #         markersize=4)


# https://stackoverflow.com/questions/25812255/row-and-column-headers-in-matplotlibs-subplots
# https://stackoverflow.com/questions/27426668/row-titles-for-matplotlib-subplot


def make_figure(df, output):
    mthds = ["treeqmc_n0_v1.0.0",
             "treeqmc_n1_v1.0.0",
             "treeqmc_n2_v1.0.0",
             "wqfm_v1.3",
             "fastral",
             "astral_3_v5.7.7"]
    mnams = [r"TREE-QMC-n0",
             r"n1",
             r"n2",
             r"wQFM",
             r"FASTRAL",
             r"ASTRAL3"]
    m = len(mthds)

    fig = plt.figure(figsize=(8, 4.75))
    gs = gridspec.GridSpec(2,3)
    ax00 = plt.subplot(gs[0,0])  ## avian 0.5X - species tree error
    ax01 = plt.subplot(gs[0,1])  ## avian 1X - species tree error
    ax02 = plt.subplot(gs[0,2])  ## avian 2X - species tree error

    ax10 = plt.subplot(gs[1,0])  ## avian 0.5X - quartet score
    ax11 = plt.subplot(gs[1,1])  ## avian 1X - quartet score
    ax12 = plt.subplot(gs[1,2])  ## avian 2X - quartet score

    rows = [ax00, ax01, ax02]
    cols = [ax10, ax11, ax12]

    scals = ["0.5X", "1X", "2X"]

    # Plot species tree error
    for i, scal in enumerate(scals):
        print(scal)

        ax = rows[i]
        varys = ["500", "true"]
        vnams = ["estimated - 500bp", "true"]
        n = len(varys)

        sers = [None] * n
        nrps = [None] * n
        for j, vary in enumerate(varys):
            sers[j] = []
            nrps[j] = []
            for k, mthd in enumerate(mthds):
                ydf = df[(df["SCAL"] == scal) &
                         (df["NBPS"] == vary) &
                         (df["MTHD"] == mthd)]
                #print(ydf)
                ydf = ydf.sort_values(by=["REPL"], ascending=True)

                #ser = ydf.SERF.values
                ser = ydf.SERF.values * 100
                sers[j].append(list(ser))
                nrps[j].append(len(ser))

        xs = []
        ys = []
        inds = varys
        for ind, ser in zip(inds, sers):
            if ser != []:
                xs.append(ind)
                ys.append(ser)
        labs = xs
        sers = ys

        xminor = []
        xmajor = []
        base = numpy.arange(1, m + 1) 
        for j, vary in enumerate(varys):
            pos = base + ((m + 1) * j)
            xminor = xminor + list(pos)
            xmajor = xmajor + [numpy.mean(pos)]

            bp = ax.boxplot(sers[j], positions=pos, widths=0.75,
                            showfliers=False, 
                            showmeans=True,
                            patch_artist=True)
            setBoxColors(bp)

        # Set labels
        #fig.suptitle("Impact of GTEE on Avian (48 taxa, 1000 gene trees)", fontsize=12)
        #ax.set_title(letters[i] + str(" Species tree height %s" % scal),
        #            loc="left", x=0.0, y=1.0, fontsize=10.5)
        ax.set_title(upperletters[i] + str(" Species tree scale %s" % scal),
                 loc="left", fontsize=11)

        if i == 0:
            ax.set_ylabel(r"Percent RF Error", fontsize=11)  

        ax.tick_params(axis='x', labelsize=10)
        ax.tick_params(axis='y', labelsize=9)
        ax.set_xlim(xminor[0]-1, xminor[-1]+1)
        ax.set_xticks(xmajor)

        test = []
        for j, vnam in enumerate(vnams):
            if nrps[j][0] != 20:
                sys.stdout.write("Avian %s - Found %d replicates!\n" % (vnam, nrps[j][0]))
            test.append(r"%s" % vnam)
        ax.set_xticklabels(test)

        if scal == "0.5X":
            yticks = [0, 5, 10, 15, 20]
        elif scal == "1X":
            yticks = [0, 3, 6, 9, 12, 15, 18]
        else:
            yticks = [0, 3, 6, 9, 12]
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


    # Plot quartet scores
    for i, scal in enumerate(scals):
        print(scal)

        ax = cols[i]
        varys = ["500", "true"]
        vnams = ["estimated - 500bp", "true"]
        n = len(varys)

        sers = [None] * n
        nrps = [None] * n
        for j, vary in enumerate(varys):
            sers[j] = []
            nrps[j] = []
            for k, mthd in enumerate(mthds):
                ydf = df[(df["SCAL"] == scal) &
                         (df["NBPS"] == vary) &
                         (df["MTHD"] == mthd)]
                ydf = ydf.sort_values(by=["REPL"], ascending=True)

                zdf = df[(df["SCAL"] == scal) &
                         (df["NBPS"] == vary) &
                         (df["MTHD"] == "TRUE")]
                zdf = zdf.sort_values(by=["REPL"], ascending=True)

                #ser = ((ydf.NQSC.values - zdf.NQSC.values) / zdf.NQSC.values) * 100
                ser = (ydf.NQSC.values - zdf.NQSC.values) * 1000
                sers[j].append(list(ser))
                nrps[j].append(len(ser))

        xs = []
        ys = []
        inds = varys
        for ind, ser in zip(inds, sers):
            if ser != []:
                xs.append(ind)
                ys.append(ser)
        labs = xs
        sers = ys

        xminor = []
        xmajor = []
        base = numpy.arange(1, m + 1) 
        for j, vary in enumerate(varys):
            pos = base + ((m + 1) * j)
            xminor = xminor + list(pos)
            xmajor = xmajor + [numpy.mean(pos)]

            bp = ax.boxplot(sers[j], positions=pos, widths=0.75,
                            showfliers=False, 
                            showmeans=True,
                            patch_artist=True)
            setBoxColors(bp)

        # Set labels
        #ax.set_title(letters[i+2] + str(" Species tree height %s" % scal),
        #            loc="left", x=0.0, y=1.0, fontsize=10.5)
        ax.set_title(upperletters[i+3] + str(" Species tree scale %s" % scal),
                 loc="left", fontsize=11)
        if i == 0:
            ax.set_ylabel(r"$\Delta$ Quartet Score $\times$ 1000", fontsize=11)  

        ax.tick_params(axis='x', labelsize=10)
        ax.tick_params(axis='y', labelsize=9)
        ax.set_xlabel("Gene Trees", fontsize=12)
        ax.set_xlim(xminor[0]-1, xminor[-1]+1)
        ax.set_xticks(xmajor)

        test = []
        for j, vnam in enumerate(vnams):
            if nrps[j][0] != 20:
                sys.stdout.write("Avian %s - Found %d replicates!\n" % (vnam, nrps[j][0]))
            test.append(r"%s" % vnam)
        ax.set_xticklabels(test)

        if scal == "0.5X":
            yticks = [-0.9, -0.6, -0.3, 0.0, 0.3, 0.6, 0.9, 1.2, 1.5]
        elif scal == "1X":
            yticks = [-0.9, -0.6, -0.3, 0.0, 0.3, 0.6, 0.9, 1.2, 1.5]
        else:
            yticks = [-0.9, -0.6, -0.3, 0.0, 0.3, 0.6, 0.9, 1.2]
        ydraw = yticks
        ymax = yticks[-1]
        ymin = -0.05 * ymax
        ax.set_ylim(ymin, ymax)
        ax.set_yticks(yticks)
    
        xs = numpy.arange(xminor[0]-1, xminor[-1]+2)
        for y in ydraw:
            if y == 0.0:
                ax.plot(xs, [y] * len(xs), "-", 
                        lw=1.5, color="black", alpha=0.3)
            else:
                ax.plot(xs, [y] * len(xs), "--", dashes=(3, 3),
                        lw=0.5, color="black", alpha=0.3)

        # Set plot axis parameters
        ax.tick_params(axis=u'both', which=u'both',length=0) # removes tiny ticks
        ax.get_xaxis().tick_bottom() 
        ax.get_yaxis().tick_left()

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    # Add legend at bottom
    gs.tight_layout(fig, rect=[0, 0.05, 1, 1])
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
              ncol=6, 
              fontsize=9.5,
              loc='lower center', 
              bbox_to_anchor=(-0.725, -0.55, 0, 3))

    # Save plot
    plt.savefig(output, format='pdf', dpi=300)

# Read and plot data
df = pandas.read_csv("../../csvs/data-all-error-and-timings.csv")
df = df[df["NGEN"] == 1000]

title = "supplot-rq2b-ils-avian-with-qscore.pdf"
make_figure(df, title)
