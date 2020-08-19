import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math
import sys

myfile = sys.argv[1] 
lengths_file=sys.argv[2]
out_file = sys.argv[3]

y_l = 6
y_u = 20
dot_norm = 15
padding=100000000
axfont = 15
threshold = 11

data = pd.read_csv(myfile, header = None, sep = "\t")

figure=None
figure = plt.figure(figsize=(12, 7),frameon=True)
ax = figure.add_subplot(111)
ax.xaxis.set_ticks_position("none")
ax.yaxis.set_ticks_position("left")
ax.set_ylabel(r'$-{log}_{10}(p)$',weight="bold")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["bottom"].set_visible(False)


ticks=[]
annots=[]
lengths_dict = {}

with open(lengths_file,'r') as f:
	count=0
	xmin = 0
	for l in f:
		count += 1
		lvals = l.strip().split()
		chrom = lvals[0]
		length = int(lvals[1])
		lengths_dict[chrom] = xmin
		xmax=xmin+length
		if (count+1) % 2 == 0:
			ax.axvspan(xmin=xmin, xmax=xmax, color="#E5E5E5")
		else:
			ax.axvspan(xmin=xmin, xmax=xmax, color="#c5c5c5")
		ticks.append((xmin + xmax) / 2)
		annots.append("chr"+str(count)+"D")
		xmin=xmax
	tot_len = xmax

data = data[data[0].isin(lengths_dict.keys())]

ax1 = ax.twinx()
ax1.xaxis.set_ticks_position("none")
ax1.yaxis.set_ticks_position("left")
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)
ax1.spines["bottom"].set_visible(False)


data['x'] = data.apply(lambda row: lengths_dict[row[0]]+row[1], axis=1)

ax1.scatter(data['x'][data[4]>=0], data[3][data[4]>=0], marker='o',c='#005AB5', s=np.ceil(data[5][data[4]>=0]/dot_norm))
ax1.scatter(data['x'][data[4]<0], data[3][data[4]<0], marker='o',c='#DC3220', s=np.ceil(data[5][data[4]<0]/dot_norm))



ax.set_xticks(ticks)
ax.set_xticklabels(annots)

ax.set_ylim([y_l,y_u])
ax1.set_ylim([y_l,y_u])

ax.set_xlim([-padding,tot_len])
ax1.set_xlim([-padding,tot_len])

for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(axfont)

for item in ([ax1.title, ax1.xaxis.label, ax1.yaxis.label] +
             ax1.get_xticklabels() + ax1.get_yticklabels()):
    item.set_fontsize(axfont)

ax.axhline(y=threshold, color='black', linestyle='--')
ax1.axhline(y=threshold, color='black', linestyle='--')

plt.show()
	
figure.savefig(out_file, format = 'eps')

