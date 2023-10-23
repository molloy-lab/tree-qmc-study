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
# Light blue - (23, 190, 207), (158, 218, 229),


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


def make_figure(df, gtre, ngen, output):
    mthds = ["treeqmc_n0_v1.0.0",
             "treeqmc_n1_v1.0.0",
             "treeqmc_n2_v1.0.0",
             "fastral",
             "astral_3_v5.7.7"]
    names = [r"TREE-QMC-n0",
             r"n1",
             r"n2",
             r"FASTRAL",
             r"ASTRAL3"]
    m = len(mthds)

    #fig = plt.figure(figsize=(8, 10))  # x,y
    fig = plt.figure(figsize=(8, 5))  # x,y
    gs = gridspec.GridSpec(2, 1)  # nrows, ncols
    ax00 = plt.subplot(gs[0, 0])  ## changing ILS
    ax10 = plt.subplot(gs[1, 0])

    # Plot Error for model condition (ILS and speciation rate)
    ax = ax00
    modls = ["0.25X\n(deep)",
             "0.25X\n(shallow)",
             "1X\n(deep)", 
             "1X\n(shallow)", 
             "5X\n(deep)",
             "5X\n(shallow)"]
    n = len(modls)

    sers = [None] * n
    nrps = [None] * n
    for j, modl in enumerate(modls):
        print(modl)

        if j < 2:
            strht = 500000 
        elif j > 3:
            strht = 10000000
        else:
            strht = 2000000

        if (j % 2) == 0:
            srate = 0.0000001  
        else:
            srate = 0.000001

        sers[j] = []
        nrps[j] = []
        for k, mthd in enumerate(mthds):
            ydf = df[(df["NTAX"] == 200) &
                     (df["STRHT"] == strht) &
                     (df["SRATE"] == srate) &
                     (df["MTHD"] == mthd)]
            #print(ydf)
            ydf = ydf.sort_values(by=["REPL"], ascending=True)

            #ser = ydf.SERF.values
            ser = ydf.SERF.values * 100
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

        bp = ax.boxplot(sers[j], positions=pos, widths=0.75,
                        showfliers=False, 
                        showmeans=True,
                        patch_artist=True)
        setBoxColors(bp)

    # Set labels
    #ax.set_title(letters[0],
    #             loc="left", x=0.0, y=1.0, fontsize=10.5)
    ax.set_title(upperletters[0],
                 loc="left", fontsize=11)
    ax.set_ylabel(r"Percent RF Error", fontsize=11)  

    ax.tick_params(axis='x', labelsize=10)
    ax.tick_params(axis='y', labelsize=9)
    #ax.set_xlabel("Species Tree Height (speciation)", fontsize=11)
    ax.set_xlim(xminor[0]-1, xminor[-1]+1)
    ax.set_xticks(xmajor)

    test = []
    for j, modl in enumerate(modls):
        if nrps[j][0] != 50:
            sys.stdout.write("%s - Found %d replicates!\n" % (modl, nrps[j][0]))
        test.append(r"%s" % modl)

    ax.set_xticklabels(test)

    if gtre == "estimated":
        if ngen == 250:
            mytitle = "Impact of ILS (200 taxa, 250 estimated gene trees)"
            yticks = list(range(0, 25+1, 5))
            ydraw = [0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25]
            ymax = 20.5
        else:
            mytitle = "Impact of ILS (200 taxa, 1000 estimated gene trees)"
            yticks = list(range(0, 18+1, 3))
            ydraw = yticks
            ymax = yticks[-1]
    elif gtre == "true":
        if ngen == 250:
            mytitle = "Impact of ILS (200 taxa, 250 true gene trees)"
            yticks = list(range(0, 15, 2))
            ydraw = yticks
            ymax = yticks[-1]
        else:
            mytitle = "Impact of ILS (200 taxa, 1000 true gene trees)"
            yticks = list(range(0, 10, 1))
            ydraw = yticks
            ymax = yticks[-1]
    
    #ax.set_title(mytitle, 
    #             loc="center", y=1.15, fontsize=12)

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

    # Plot runtime for varying ILS
    ax = ax10
    #ax.set_title(letters[1],
    #             loc="left", x=0.0, y=1.0, fontsize=10.5)
    ax.set_title(upperletters[1],
                 loc="left", fontsize=11)


    for k, mthd in enumerate(mthds):
        av = []
        se = []
        for j, modl in enumerate(modls):
            print(modl)

            if j < 2:
                strht = 500000
            elif j > 3:
                strht = 10000000
            else:
                strht = 2000000

            if (j % 2) == 0:
                srate = 0.0000001 
            else:
                srate = 0.000001
        
            xdf = df[(df["NTAX"] == 200) &
                     (df["STRHT"] == strht) &
                     (df["SRATE"] == srate) &
                     (df["MTHD"] == mthd)]

            if xdf.shape[0] > 0:
                vals = xdf.SECS.values / (60.0 * 60.0)

                av.append(numpy.mean(vals))
                se.append(numpy.std(vals) / numpy.sqrt(vals.size))
    
        av = numpy.array(av)
        se = numpy.array(se)

        labls = [1, 2, 3, 4, 5, 6] #modls
        ax.plot(labls, av, '-', color=tableau20[2*k], lw=2.5, 
                    label=names[k])
        ax.fill_between(labls, av - se, av + se, 
                        color=tableau20[2*k+1], alpha=0.5)

    ax.set_xticks([1, 2, 3, 4, 5, 6])
    ax.set_xticklabels(modls)

    if gtre == "estimated":
        if ngen == 250:
            yticks = [0.0, 0.05, 0.1, 0.15, 0.2]
            ydraw = yticks
            ymax = yticks[-1]
        else:
            yticks = [0.25, 0.5, 0.75, 1.0, 1.25, 1.5]
            ydraw = yticks
            ymax = yticks[-1]
    elif gtre == "true":
        if ngen == 250:
            yticks = [0.0, 0.02, 0.04, 0.06, 0.08, 0.1]
            ydraw = yticks
            ymax = yticks[-1]
        else:
            yticks = [0.0, 0.2, 0.4, 0.6, 0.8, 1]
            ydraw = yticks
            ymax = yticks[-1]

    ax.set_ylim(0, ymax)
    ax.set_yticks(yticks)
    xs = ax.xaxis.get_majorticklocs()
    for y in ydraw:        
        ax.plot(xs, [y] * len(xs), "--", dashes=(3, 3),
                lw=0.5, color="black", alpha=0.3)

    ax.set_xlabel(r"Species Tree Height (speciation)", fontsize=12)
    ax.set_ylabel(r'Runtime (h)', fontsize=11)
    ax.tick_params(axis='x', labelsize=10)
    ax.tick_params(axis='y', labelsize=9)

    ax.get_xaxis().tick_bottom() 
    ax.get_yaxis().tick_left() 
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


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
                  frameon=False, ncol=5, 
                  fontsize=9.5,
                  loc='upper center', 
                  bbox_to_anchor=(0.5, -0.465))

    # Save plot
    #gs.tight_layout(fig, rect=[0, 0, 1, 1])
    plt.savefig(output, format='pdf', dpi=300)

# Read and plot data
for i in range(4):
    df = pandas.read_csv("../../csvs/data-all-error-and-timings-ex.csv")

    if i < 2:
        gtre = "estimated"
        df = df[df["GTRE"] == gtre + "genetre"]

        if i == 0:
            tag = "supplot"
            ngen = 250
        else:
            tag = "plot"
            ngen = 1000
    else:
        gtre = "true"
        df = df[df["GTRE"] == gtre + "genetrees"]

        if i == 2:
            tag = "supplot"
            ngen = 250
        else:
            tag = "supplot"
            ngen = 1000

    df = df[df["NGEN"] == ngen]

    title = str("%s-rq2a-ils-%d-%sgtres.pdf" % (tag, ngen, gtre))
    make_figure(df, gtre, ngen, title)
