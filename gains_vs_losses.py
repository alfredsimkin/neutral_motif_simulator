import custom
import sys
data_folder=sys.argv[1]

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
	plt.plot((0, 65536), (3277, 3277), 'b', linestyle='dashed')
	plt.plot((3277, 3277), (0, 65536), 'b', linestyle='dashed')
	m, b, r, linear_p_value, std_err=scipy.stats.linregress(first, second)
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
#	plt.xlabel(first_label, size=15)
#	plt.ylabel(second_label, size=15)
#	linear_p_value=numpy.round(linear_p_value, )
	plt.title('m='+str(custom.round_scientific(m))+' b='+str(custom.round_scientific(b))+' r2='+str(custom.round_scientific(r**2))+' p='+str(custom.round_scientific(linear_p_value)), size=15)
#	plt.grid(True, linestyle='-')
	ax = plt.gca()
	ax.tick_params(direction='out')
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.yaxis.set_ticks_position('left')
	ax.xaxis.set_ticks_position('bottom')
	ax.xaxis.set_tick_params(width=2)
	ax.yaxis.set_tick_params(width=2)
	plt.savefig(data_folder+first_label+'_vs_'+second_label+'.pdf', bbox_inches='tight')
	plt.clf()

gains=sorted([line.strip().split('\t') for line in open(data_folder+'mir_ranks_gain')])
losses=sorted([line.strip().split('\t') for line in open(data_folder+'mir_ranks_loss')])
output_file=open(data_folder+'gains_vs_losses_raw_data', 'w')
output_file.write('8mer\tmir_name\tgain_rank\tloss_rank\n')
for line_number, line in enumerate(gains):
	output_file.write('\t'.join([line[0], line[5], line[1], losses[line_number][1]])+'\n')
gains=[float(line[1]) for line in gains]
losses=[float(line[1]) for line in losses]
scatter_it(gains, losses, 'miRNA_gain_full_rank', 'miRNA_loss_full_rank_test', lims=[[-50,66000], [-50,66000]], faded=False)
#gain_ranks, loss_ranks=custom.make_ranks(gains, losses)
#custom.scatter_it(gain_ranks, loss_ranks, 'miRNA_gain_rel_rank', 'miRNA_loss_rel_rank_test', lims=False, faded=False)
