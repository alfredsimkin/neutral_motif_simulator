'''
Strings together several scripts to automate analysis. User inputs folder where
simulated counts are located, replicate start, replicate end to be analyzed,
folder where stats and graphs should go (Doesn't have to exist), input file of
real counts. sig_miRNAs.py uses a default miRNA file that can be swapped out,
and graph_distributions_matplotlib.py uses a default 'normal distribution' that
will also need to be swapped out if the x-axis scale of the graph changes.
This version handles all possible seed conversion events
'''

import subprocess

#still need to be integrated

data_folder='/mnt/hulk_tank/alfred_stuff/silver_surfer/Grimson_new_phylogeny/data/other_data/'
rep_start, rep_end='0','100'
#program expects seeds to be listed from smallest to largest
seed_types='none_6mer_7mer1A_7merm8_8mer'
miRNA_file='/mnt/hulk_tank/alfred_stuff/silver_surfer/Grimson_new_phylogeny/data/real_eightmers_human_mouse_fish_best_group3'
stats_folder=data_folder+'_'.join(['conversions', rep_start, rep_end, seed_types, 'human_mouse_fish_group3/'])
subprocess.call(['mkdir', stats_folder])
import subprocess
seed_comparisons=[]


#print 'counting real events'
subprocess.call(['python', 'count_every_real_thing.py', data_folder, stats_folder, seed_types])

#print 'counting simulated events'
subprocess.call(['python', 'count_every_simulated_thing.py', data_folder, seed_types, rep_start, rep_end])

seed_types=seed_types.split('_')
for small_number, small_seed in enumerate(seed_types):
	for large_seed in seed_types[small_number+1:]:
		seed_comparisons.append(small_seed+'_'+large_seed)

for seed_comparison in seed_comparisons:
	print '*********\nseed is', seed_comparison, '\n*********'
	print 'making distributions'
	subprocess.call(['python', 'make_every_distribution.py', stats_folder, rep_start, rep_end, seed_comparison])
#ranks miRNA gain and loss rates relative to other 8mers and send to files. Not
#used by other programs in this collection, more useful to other pipelines
	print 'ranking miRNAs'
	subprocess.call(['python', 'extract_every_miRNA_rank.py', stats_folder, seed_comparison, miRNA_file])

#plots the rates from extract_miRNA_ranks.py
	print 'plotting miRNA gains against losses'
	subprocess.call(['python', 'every_gains_vs_losses.py', stats_folder, seed_comparison])
	
#bins all 8mers and conserved miRNAs by their observed number of std devs away
#from simulated values, calculates stats for each bin
	print 'binning all 8mers and miRNA 8mers into stdev bins'
	subprocess.call(['python', 'bin_every_distribution.py', stats_folder, seed_comparison])
	
	print 'graphing the bins'
	subprocess.call(['python', 'graph_every_distribution_matplotlib.py', stats_folder, seed_comparison])
