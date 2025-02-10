'''
input files: 
list of 8mers of interest, conservation levels
output files: conservation_vs_score_(gains|losses|sum)

this version only plots miRNAs conserved above a fixed threshold.
'''
import custom
from decimal import Decimal

threshold=1
data_folder='/mnt/hulk_tank/alfred_stuff/silver_surfer/Grimson_new_phylogeny/data/other_data/paper_stats/'
rank_file=data_folder+'mir_ranks_'
mir_file='/mnt/hulk_tank/alfred_stuff/silver_surfer/Grimson_new_phylogeny/data/real_eightmers_human_mouse_fish_best'

super_mirs={line.strip().split('\t')[0]:line.strip().split('\t')[1] for line in open(mir_file)}

stat_list=['gain', 'loss']
cons_dict={}

def scatter_it(first, second, first_label, second_label, lims=False, faded=True):
	'''
	takes a list of x values ('first') and y values ('second'), along with x
	labels and y labels, and plots them as a scatter plot. Unless specified,
	limits are auto generated, and data points are faded. Figure will be saved
	as x_label_vs_y_label.png. A linear regression line is also plotted and
	listed as the title of the graph.
	'''
	import scipy.stats
	import matplotlib as mpl
	mpl.use('Agg') #comment this in if using a nongraphical interface, comment out to enable graphical interface
	import matplotlib.pyplot as plt
	import numpy
	m, b, r, linear_p_value, std_err=scipy.stats.linregress(first, second)
#	r, rank_p_value=scipy.stats.spearmanr(first, second)
	y=[m*x+b for x in first]
	if lims!=False:
		plt.xlim(lims[0])
		plt.ylim(lims[1])
	if faded:
		plt.scatter(first, second, color='red', linewidth=0, alpha=0.08)
	else:
		plt.scatter(first, second, color='red', linewidth=0)		
	plt.plot(first, y, color="black")
#	plt.tick_params(labelsize=17)
	plt.xlabel(first_label)
	plt.ylabel(second_label)
	ax = plt.gca()
	ax.xaxis.set_tick_params(width=2)
	ax.yaxis.set_tick_params(width=2)
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.yaxis.set_ticks_position('left')
	ax.xaxis.set_ticks_position('bottom')
	ax.tick_params(direction='out')
#	plt.xlabel(first_label, size=15)
#	plt.ylabel(second_label, size=15)
#	linear_p_value=numpy.round(linear_p_value, )
	plt.title('m='+str(custom.round_scientific(m))+' b='+str(custom.round_scientific(b))+' r2='+str(custom.round_scientific(r**2))+' p='+str(custom.round_scientific(linear_p_value)), size=15)
#	plt.title('m='+str(custom.round_scientific(m))+' b='+str(custom.round_scientific(b))+' r2='+str(custom.round_scientific(r**2))+' p='+str(custom.round_scientific(rank_p_value)), size=15)
#	plt.grid(True, linestyle='-')
	plt.savefig(first_label+'_vs_'+second_label+'.pdf', bbox_inches='tight')
	plt.clf()


def make_ranks(x, y):
	import scipy
	x=scipy.stats.stats.rankdata(x)
	y=scipy.stats.stats.rankdata(y)
	return x, y

for line in open('mature_mir_branch_lengths'):
	line=line.strip().split()
	if line[1] not in cons_dict or float(line[2])>cons_dict[line[1]]:
		cons_dict[line[1]]=float(line[2])

for stat in stat_list:
#	super_ranks={line.strip().split('\t')[0]:int(line.strip().split('\t')[-1]) for line in open(rank_file+stat)}
	super_ranks={line.strip().split('\t')[0]:int(line.strip().split('\t')[1]) for line in open(rank_file+stat)}
	x_values, y_values, super_x_values, super_y_values, super_z_values=[],[],[],[],[]
	output_file=open('conservation_vs_deep_stdev_'+stat, 'w')
	output_file.write('\t'.join(['8mer', 'std_devs', 'cons_level', 'status'])+'\n')
	for line_number, line in enumerate(open(data_folder+'/real_'+stat)):
		line=line.strip().split('\t')
#		eightmer, score=line[0], line_number
		eightmer, score=line[0], float(line[1])
		if eightmer in cons_dict:
			if eightmer in super_mirs:
				status=super_mirs[eightmer]
				super_x_values.append(cons_dict[eightmer])
				super_y_values.append(score)
				super_z_values.append(super_ranks[eightmer])
			else:
				status='not_super'
			if cons_dict[eightmer]>threshold:
				x_values.append(cons_dict[eightmer])
				y_values.append(score)
				string_list=map(str, [eightmer, score, cons_dict[eightmer], status])
				output_file.write('\t'.join(string_list)+'\n')
	scatter_it(x_values, y_values, str(stat)+'_conservation_level', 'deep_stdev_value2_final', [[0,10],[-37, 12]], False)
	output_file.close()
