'''
Strings together several scripts to automate analysis.

This version should be run after running all_together.py. It redoes the
analysis portion of all_together.py, calling a script that locates sites that
are not perfectly conserved and only considering 8mers whose nucleotides
overlap these sites
'''


import subprocess
import tree_manips
import time

data_folder='/mnt/hulk_tank/alfred_stuff/silver_surfer/Grimson_new_phylogeny/data/other_data/'
miRNA_file='/mnt/hulk_tank/alfred_stuff/silver_surfer/Grimson_new_phylogeny/data/real_eightmers_human_mouse_fish_best'
#maf_folder='/media/alfred/data_drive/big_data/reference_files/multiz_100way/'
#target_species=['human', 'ponAbe2', 'saiBol1', 'speTri2', 'turTru2', 'dasNov3', 'loxAfr3', 'susScr3', 'cavPor3', 'rn5']
#common_sp=['human', 'orangutan', 'squirrel monkey', 'squirrel', 'dolphin', 'armadillo', 'elephant', 'pig', 'guinea pig', 'rat']
#input_bed_file='4-28-15_refseq.bed'
#input_tree='hg19.100way.nh'
#bed_prefix=input_bed_file[:-4]
rep_start, rep_end='0','100'

stats_folder=data_folder+'paper_test_perfect_8mers/'
subprocess.call(['mkdir', stats_folder])
mask_file=stats_folder+'single_nucleotide_substitutions'

print 'making single nucleotide substitution mask'
subprocess.call(['python', 'locate_SNSs.py', data_folder, mask_file])

print 'counting real turnover events'
subprocess.call(['python', 'count_real_minus_perfect.py', data_folder, stats_folder, mask_file])

#calculates every real 8mer's standard deviations relative to the mean of its
#simulated values, outputs to stats folder. Very slow, memory intensive
print 'making distributions'
subprocess.call(['python', 'make_distributions_minus_perfect_expanded.py', stats_folder, rep_start, rep_end, mask_file])

print 'pulling out the std devs for conserved miRNAs relative to the full dataset'
subprocess.call(['python', 'sig_miRNAs.py', stats_folder, miRNA_file])

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
