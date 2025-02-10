import math
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy
import sys

stats_folder, seed_comparison=sys.argv[1:]
#data_folder='/media/alfred/data_drive/big_data/turnover_rates/no_windows_revised/stats/first_100_replicates/'

#normal_list=[float(line.strip().split()[1]) for line in open('/media/alfred/data_drive/big_data/turnover_rates/normal_dist')]
stat_list=['gain', 'loss', 'sum']
window_size=2

def window_average(window_size, input_x, input_y):
	'''
	returns a window average of y_values for each x point that is at least window_size positions from an edge
	'''
	new_x, new_y=[],[]
	for x_number, x_value in enumerate(input_x):
		if x_number>=window_size and x_number<len(input_x)-window_size:
			new_x.append(x_value)
			to_average=0
			for value in range(x_number-window_size, x_number+window_size+1):
				to_average+=input_y[value]
			new_y.append(to_average/(2*window_size+1))
	return new_x, new_y

def print_smoothed(x_smoothed, mir_smoothed, background_smoothed, stat):
	smoothed_output=open(stats_folder+'distribution_'+stat+'_'+seed_comparison+'_smoothed', 'w')
	smoothed_output.write('std_devs\tmir_smoothed\tfull_smoothed\n')
	for x_number, x in enumerate(x_smoothed):
		smoothed_output.write(str(x)+'\t'+str(mir_smoothed[x_number])+'\t'+str(background_smoothed[x_number])+'\n')

for stat in stat_list:
	print 'plotting', stat
	input_list=[line.strip().split() for line in open(stats_folder+'distribution_'+stat+'_'+seed_comparison)]
	std_devs, real, full, colors=[],[],[],[]
	for line in input_list[1:]:
		std_devs.append(float(line[0]))
		real.append(float(line[3]))
		full.append(float(line[4]))
		colors.append(line[6])
	new_std_devs, new_full=window_average(window_size, std_devs, full)
	plt.plot(new_std_devs, new_full, c="blue")
	new_std_devs, new_real=window_average(window_size, std_devs, real)
	print_smoothed(new_std_devs, new_real, new_full, stat)
	plt.plot(new_std_devs, new_real, c="red")
#	plt.plot(std_devs, full, marker='.')
#	plt.plot(std_devs, normal_list)
#	plt.plot(std_devs, real, linestyle='--')
#	plt.scatter(std_devs, real, linewidth=0, c=colors)
	lims = plt.xlim([-30, 30])
	lims=plt.ylim([0,0.20])
	ax = plt.gca()
	ax.xaxis.set_tick_params(width=2)
	ax.yaxis.set_tick_params(width=2)
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.yaxis.set_ticks_position('left')
	ax.xaxis.set_ticks_position('bottom')
	ax.tick_params(direction='out')
	plt.xlabel('standard deviations', size=20)
	plt.ylabel('proportion', size=20)
#	plt.grid(True)
	plt.title(seed_comparison+'_'+stat)
#	miRNA=plt.Line2D((0,0),(0,0), color='red', linestyle='--')
#	normal=plt.Line2D((0,0),(0,0), color='green')
	all_8mers=plt.Line2D((0,0),(0,0), color='blue')
#	sig_05=plt.Line2D((0,0),(0,0), color='DarkViolet', marker='o', linestyle='')
#	sig_01=plt.Line2D((0,0),(0,0), color='orange', marker='o', linestyle='')
#	sig_001=plt.Line2D((0,0),(0,0), color='red', marker='o', linestyle='')
#	plt.legend([miRNA, normal, all_8mers, sig_05, sig_01, sig_001], ["cons. miRNAs", "normal_dist", 'all_8mers', 'sig 0.05', 'sig 0.01', 'sig 0.001'], prop={'size':8})
	#plt.legend([normal, all_8mers, miRNA], ['normal_dist', 'all_8mers', 'cons. miRNAs'], prop={'size':8})
	plt.savefig(stats_folder+stat+'_'+seed_comparison+'_distributions4.pdf', bbox_inches='tight')
	#plt.show()
	plt.clf()
