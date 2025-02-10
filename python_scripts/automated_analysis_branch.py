'''
Strings together several scripts to automate analysis.

This version should be run after running all_together.py. It redoes the
analysis portion of all_together.py, calling scripts that count real and
simulated counts in a branch-specific manner, and then graph the relationship
between branch length and total counts
'''
import cPickle
import subprocess
import tree_manips
import time

data_folder='/mnt/hulk_tank/alfred_stuff/silver_surfer/Grimson_new_phylogeny/data/other_data/'
miRNA_file='/mnt/hulk_tank/alfred_stuff/silver_surfer/Grimson_new_phylogeny/data/real_eightmers_human_mouse_fish_best'
rep_start, rep_end='0','100'

stats_folder=data_folder+'branch_specific_stats/'
subprocess.call(['mkdir', stats_folder])

gene_dict=cPickle.load(open(data_folder+'reformatted_anc_dictU', 'rb'))
gene=gene_dict.keys()[0]
'''
print 'counting real turnover events'
subprocess.call(['python', 'count_real_branch.py', data_folder+'reformatted_anc_dictU', stats_folder, 'real_branch_counts'])

print 'counting simulated turnover events'
for rep in range(int(rep_start), int(rep_end)):
	print 'analyzing replicate', rep
	subprocess.call(['python', 'count_real_branch.py', data_folder+'simulations/formatted_sim_wrong_anc_dict'+str(rep), stats_folder, 'sim_branch_counts'+str(rep)])

#calculates every real 8mer's standard deviations relative to the mean of its
#simulated values, outputs to stats folder. Very slow, memory intensive
print 'making distributions'
subprocess.call(['python', 'make_distributions_branch_expanded.py', stats_folder, rep_start, rep_end])

#ranks miRNA gain and loss rates relative to other 8mers and send to files
print 'ranking miRNAs'
for ancestor in gene_dict[gene]:
	for descendant in gene_dict[gene][ancestor]:
		subprocess.call(['python', 'extract_miRNA_branch_ranks.py', stats_folder, miRNA_file, ancestor+'-'+descendant])
'''
#bin_distributions.py

#graph_distributions.py

#branch_length_vs_site_depletion.py

#branch_length_vs_rank_skew.py

#plots the ranks from extract_miRNA_ranks.py


#bins all 8mers and conserved miRNAs by their observed number of std devs away
#from simulated values, calculates stats for each bin
print 'binning all 8mers and miRNA 8mers into stdev bins'
for ancestor in gene_dict[gene]:
	for descendant in gene_dict[gene][ancestor]:
		subprocess.call(['python', 'bin_branch_distributions.py', stats_folder, ancestor+'-'+descendant])



'''

print 'graphing the bins'
subprocess.call(['python', 'graph_distributions_matplotlib.py', stats_folder])
'''

#this script is currently 'stand alone' and should be integrated into this automated analysis pipeline:
#branch_length_vs_stats.py
