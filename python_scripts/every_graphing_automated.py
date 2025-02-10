'''
Strings together several scripts to automate analysis. User inputs folder where
simulated counts are located, replicate start, replicate end to be analyzed,
folder where stats and graphs should go (Doesn't have to exist), input file of
real counts. sig_miRNAs.py uses a default miRNA file that can be swapped out,
and graph_distributions_matplotlib.py uses a default 'normal distribution' that
will also need to be swapped out if the x-axis scale of the graph changes.
This version handles all possible seed conversion events
'''
input_folder='/media/alfred/data_drive/big_data/turnover_rates/4-28-15_deeper_phylogeny/'
rep_start, rep_end='0','100'
output_folder=input_folder+'stats/'
import subprocess
seed_types=['none', '6mer', '7mer1A', '7merm8', '8mer']
seed_comparisons=[]
for small_number, small_seed in enumerate(seed_types):
	for large_seed in seed_types[small_number+1:]:
		seed_comparisons.append(small_seed+'_'+large_seed)
for seed_comparison in seed_comparisons:
	print '*********\nseed is', seed_comparison, '\n*********'
	print 'making distributions'
	subprocess.call(['python', 'make_every_distribution.py', input_folder, rep_start, rep_end, seed_comparison])

	print 'pulling out the std devs for conserved miRNAs relative to the full dataset'
	subprocess.call(['python', 'sig_miRNAs.py', output_folder])

#ranks miRNA gain and loss rates relative to other 8mers and send to files. Not
#used by other programs in this collection, more useful to other pipelines
	print 'ranking miRNAs'
	subprocess.call(['python', 'extract_every_miRNA_rank.py', output_folder, seed_comparison])

#plots the rates from extract_miRNA_ranks.py
	print 'plotting miRNA ranks'
	subprocess.call(['python', 'gains_vs_losses.py', output_folder]) #alexa note: this is dispensable, it should produce graphs similar to mine, but you'll probably want to make your own graphing script

#bins all 8mers and conserved miRNAs by their observed number of std devs away
#from simulated values, calculates stats for each bin
	print 'binning all 8mers and miRNA 8mers into stdev bins'
	subprocess.call(['python', 'bin_every_distribution.py', output_folder, seed_comparison])

	print 'graphing the bins'
	subprocess.call(['python', 'graph_every_distribution_matplotlib.py', output_folder, seed_comparison])
