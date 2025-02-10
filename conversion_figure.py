'''
Converts a table of p-values into a heatmap figure. Thanks to Katla for
thinking of this.
'''
import math
state, u, p, size, one, two, direction=range(7)
labels=['none', '6mer', '7mer1A', '7merm8', '8mer']

input_data=[[0]*len(labels) for number in range(len(labels))]
for line in open('mannu_conversions'):
	line=line.strip().split('\t')
	if line[direction]=='gain':
		x, y=line[one], line[two]
	else:
		x, y=line[two], line[one]
	p_value=float(line[p])
	if state=='x_smaller':
		p_value=1-p_value
		line[u]=float(line[size])-float(line[u])
	p_value=-1*math.log(p_value, 10)
	skew=1-float(line[u])/float(line[size])
	x, y=labels.index(x), labels.index(y)
	input_data[y][x]=p_value
	#input_data[y][x]=skew


import matplotlib as mpl
mpl.use('Agg') #comment this in if using a nongraphical interface, comment out to enable graphical interface
import matplotlib.pyplot as plt
import numpy as np
column_labels = labels
row_labels = labels
data=np.array(input_data)
fig, ax = plt.subplots()
#heatmap = ax.pcolor(data, vmin=0.5, vmax=1.0, cmap=plt.cm.gist_earth_r)
heatmap = ax.pcolor(data, vmin=1.3, vmax=36, cmap=plt.cm.gist_earth_r)
plt.colorbar(heatmap)
# put the major ticks at the middle of each cell
ax.set_xticks(np.arange(data.shape[0])+0.5, minor=False)
ax.set_yticks(np.arange(data.shape[1])+0.5, minor=False)

# want a more natural, table-like display
ax.invert_yaxis()
ax.xaxis.tick_top()

xlabels = ax.get_xticklabels()
ylabels = ax.get_yticklabels()
plt.setp(xlabels, rotation=45, fontsize=15)
plt.setp(ylabels, fontsize=15)

ax.set_xticklabels(row_labels, minor=False)
ax.set_yticklabels(column_labels, minor=False)
ax.tick_params(axis=u'both', which=u'both',length=0)
#plt.figure(figsize=(20,20))
#plt.subplots_adjust(left=0.2, right=0.9, top=0.8, bottom=0.1)
plt.savefig('all_p_test_all.pdf', bbox_inches='tight')
