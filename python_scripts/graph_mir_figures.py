'''
Smooths points and graphs the output of populate_mir_figures.py
'''

import cPickle
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy
from matplotlib import colors
import sys

figure_folder, smoothing_factor=sys.argv[1:]
populated_figures=cPickle.load(open(figure_folder+'populated_figures', 'rb'))
k=float(smoothing_factor)

def calculate_distance(first_point, second_point):
	x_distance=abs(first_point[0]-second_point[0])
	y_distance=abs(first_point[1]-second_point[1])
	return (x_distance**2+y_distance**2)**0.5

def distribute_values(current_x, current_y, intensity, value_list):
	for x in range(84):
		for y in range(84):
			r=calculate_distance([current_x, current_y], [x, y])
			result=(1/(r**2+(k*84)))*intensity
#			if value_list[x][y]!=0:
#				print '\n', value_list[x][y]
			value_list[x][y]+=result
#			print value_list[x][y]
	return value_list

figures=[('CGGTACGA', 'hsa-miR-126-3p')]
for figure in populated_figures:
	final_values={}
#for figure in figures:
	print figure
	tissue=populated_figures[figure][0]
	for gene_number, gene in enumerate(populated_figures[figure][1]):
		if gene_number%1000==0:
			print figure, gene_number
		for event in populated_figures[figure][1][gene]:
			if event not in final_values:
				final_values[event]=[[0]*84 for number in range(84)]
			x, y, intensity=populated_figures[figure][1][gene][event]
			#if intensity!=0:
				#print event, gene, intensity
			final_values[event]=distribute_values(x, y, intensity, final_values[event])
	for event in final_values:
		max_smoothed, min_smoothed=0,0
		for x_number, x in enumerate(final_values[event]):
			for y_number, y in enumerate(x):
				smoothed=y
				if smoothed>max_smoothed:
					max_smoothed=smoothed
				if smoothed<min_smoothed:
					min_smoothed=smoothed
		print figure, tissue, event, max_smoothed, min_smoothed
		cmap = colors.ListedColormap([(0.7, 0.1, 0.1, 0.8), (0.7, 0.1, 0.1, 0.6), (0.7, 0.1, 0.1, 0.4), (0.7, 0.1, 0.1, 0.2), (0.7, 0.1, 0.1, 0.1), (0.1, 0.1, 0.7, 0.1), (0.1, 0.1, 0.7, 0.2), (0.1, 0.1, 0.7, 0.4), (0.1, 0.1, 0.7, 0.6), (0.1, 0.1, 0.7, 0.8)])
		smallest=min(abs(max_smoothed), abs(min_smoothed))
		increment=(smallest)/5.0
		bounds=[(-1*smallest)+(value_number*increment) for value_number in range(11)]
		print bounds
		#bounds=[-0.05, -0.04, -0.03, -0.02, -0.01, 0, 0.01, 0.02, 0.03, 0.04, 0.05]
		norm = colors.BoundaryNorm(bounds, cmap.N)
		heatmap = plt.pcolor(numpy.array(final_values[event]), cmap=cmap, norm=norm)
		cPickle.dump(final_values[event], open(figure_folder+figure[1]+'_'+tissue+'_'+event+'_'+str(k)+'.pickle', 'wb'), -1)
		plt.colorbar(heatmap, ticks=bounds)
		plt.savefig(figure_folder+figure[1]+'_'+tissue+'_'+event+'_'+str(k)+'.png', bbox_inches='tight')
		plt.clf()
'''
	heatmap=plt.pcolor(numpy.array(final_values[event]))
	cbar = plt.colorbar(heatmap)
#	cbar.ax.set_yticklabels(['0','1','2','>3'])
	cbar.set_label('some_name', rotation=270)
	plt.savefig(event+'something.png', bbox_inches='tight')
	plt.clf()
'''

