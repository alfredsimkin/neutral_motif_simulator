'''
Strings together several scripts to automate analysis.

This version should be run after running all_together.py. It redoes the
analysis portion of all_together.py, calling a script that locates sites in
contexts with > or < 50% AU composition. Check for > by typing 'high' as the
high_low variable, and < by typing 'low' as the high_low variable
'''


import subprocess
import tree_manips
import time

data_folder='/mnt/hulk_tank/alfred_stuff/silver_surfer/Grimson_new_phylogeny/data/other_data/'
miRNA_file='/mnt/hulk_tank/alfred_stuff/silver_surfer/Grimson_new_phylogeny/data/real_eightmers_human_mouse_fish_best'
rep_start, rep_end='0','100'

high_low='high'
stats_folder=data_folder+high_low+'_AU_context_paper/'
subprocess.call(['mkdir', stats_folder])
mask_file=stats_folder+high_low+'_AU_context_mask'

print 'making AU mask'
subprocess.call(['python', 'locate_AU.py', data_folder, high_low, mask_file])

print 'counting real turnover events'
subprocess.call(['python', 'count_real_AU_context.py', data_folder, stats_folder, mask_file, high_low])

#calculates every real 8mer's standard deviations relative to the mean of its
#simulated values, outputs to stats folder. Very slow, memory intensive
print 'making distributions'
subprocess.call(['python', 'make_distributions_mask_AU_expanded.py', stats_folder, rep_start, rep_end, mask_file, high_low])

#ranks miRNA gain and loss rates relative to other 8mers and send to files
print 'ranking miRNAs'
subprocess.call(['python', 'extract_miRNA_ranks.py', stats_folder, miRNA_file])

#plots the ranks from extract_miRNA_ranks.py
print 'plotting miRNA ranks'
subprocess.call(['python', 'gains_vs_losses.py', stats_folder])


#bins all 8mers and conserved miRNAs by their observed number of std devs away
#from simulated values, calculates stats for each bin
print 'binning all 8mers and miRNA 8mers into stdev bins'
subprocess.call(['python', 'bin_distributions.py', stats_folder])

print 'graphing the bins'
subprocess.call(['python', 'graph_distributions_matplotlib.py', stats_folder])
