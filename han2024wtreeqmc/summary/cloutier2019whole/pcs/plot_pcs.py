import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import rc
import numpy
import pandas
import sys

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{helvet} \usepackage{sfmath}')

w10 = 10

# Read data frame
df_h = pandas.read_csv("wtreeqmc_hybrid_n2_treeA-for-pcs.tsv", sep="\t")
df_n = pandas.read_csv("wtreeqmc_none_n2_treeA-for-pcs.tsv", sep="\t")

# Normalize quartets
pos = df_h["position"].values

# Plot figure
fig = plt.figure(figsize=(8, 5))  # x,y
gs = gridspec.GridSpec(2, 1)  # nrows, ncols
ax0 = plt.subplot(gs[0, 0])
ax1 = plt.subplot(gs[1, 0])

f1_n = df_n["f_xy|zw"].values
f2_n = df_n["f_xz|yw"].values
f3_n = df_n["f_xw|yz"].values
en_n = df_n["totalf"].values
q1_n = f1_n / en_n
q2_n = f2_n / en_n
q3_n = f3_n / en_n


inds = numpy.array([x for x in range(0, len(pos))])

# Plot normal weights
windows = [x for x in range(pos[0], pos[-1], w10)]
q1dots = []
q2dots = []
q3dots = []
wdots = []
for i, win in enumerate(windows[:-1]):
	start = windows[i]
	end = windows[i+1]
	print("window %d to %d" % (start, end))

	# Now find positions between window
	find = inds[(pos >= start) & (pos < end)]
	if (len(find) > 0):
		q1dots.append(numpy.mean(q1_n[find]))
		q2dots.append(numpy.mean(q2_n[find]))
		q3dots.append(numpy.mean(q3_n[find]))
	else:
		print("    Nothing found in window")
		q1dots.append(0)
		q2dots.append(0)
		q3dots.append(0)

	wdots.append(i * w10)



ax0.set_title(r"\textbf{a)} Unweighted", loc="left", fontsize=12)

ax0.plot(wdots, q1dots, '.', color='red', alpha=0.8, markeredgewidth=0.0,
        label=r"T, K+C+E $|$ R, O+C")

ax0.plot(wdots, q2dots, '.', color='orange', alpha=0.8, markeredgewidth=0.0,
        label=r"T, R $|$ K+C+E, O+C")

ax0.plot(wdots, q3dots, '.', color='blue', alpha=0.8, markeredgewidth=0.0,
        label=r"T, O+C $|$ K+C+E, R")

ax0.axvline(x=(3158 - 105 + 1), color='k',
           label='UCEs with homology errors are on right of line')

ax0.spines["top"].set_visible(False)
ax0.spines["right"].set_visible(False)

ax0.set_ylabel("Normalized QQS")
ax0.set_ylim(0.0, 1.0)
ax0.legend(frameon=False, loc="upper left", fontsize=9, ncol=2)


# Plot hybrid
f1_h = df_h["f_xy|zw"].values
f2_h = df_h["f_xz|yw"].values
f3_h = df_h["f_xw|yz"].values
en_h = df_h["totalf"].values
q1_h = f1_h 
q2_h = f2_h 
q3_h = f3_h 


windows = [x for x in range(pos[0], pos[-1], w10)]
q1dots = []
q2dots = []
q3dots = []
wdots = []
for i, win in enumerate(windows[:-1]):
	start = windows[i]
	end = windows[i+1]
	print("window %d to %d" % (start, end))

	# Now find positions between window
	find = inds[(pos >= start) & (pos < end)]
	if (len(find) > 0):
		q1dots.append(numpy.mean(q1_h[find]))
		q2dots.append(numpy.mean(q2_h[find]))
		q3dots.append(numpy.mean(q3_h[find]))
	else:
		print("    Nothing found in window")
		q1dots.append(0)
		q2dots.append(0)
		q3dots.append(0)

	wdots.append(i * w10)

ax1.set_title(r"\textbf{b)} Weighted (hybrid)", loc="left", fontsize=12)

ax1.plot(wdots, q1dots, '.', color='red', alpha=0.8, markeredgewidth=0.0)
ax1.plot(wdots, q2dots, '.', color='orange', alpha=0.8, markeredgewidth=0.0)
ax1.plot(wdots, q3dots, '.', color='blue', alpha=0.8, markeredgewidth=0.0)
ax1.axvline(x=(3158 - 105 + 1), color='k')

ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)

ax1.set_xlabel("3158 UCE gene trees")
ax1.set_ylabel("Normalized QQS")
ax1.set_ylim(0.0, 1.0)

gs.tight_layout(fig, rect=[0, 0, 1, 1])
plt.savefig("pcs.pdf", format='pdf', dpi=300)

plt.show()
