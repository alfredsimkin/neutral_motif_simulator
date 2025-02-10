'''
Driver script to analyze the effects each gene has on overall turnover rates
'''

import subprocess

master_folder='/mnt/hulk_tank/alfred_stuff/silver_surfer/Grimson_new_phylogeny/data/other_data/'
simulation_folder=master_folder+'simulations/'
output_folder=master_folder+'single_gene_removal/best_mirs/vs_reps_test/'
miRNA_file='/mnt/hulk_tank/alfred_stuff/silver_surfer/Grimson_new_phylogeny/data/real_eightmers_human_mouse_fish_best'
subprocess.call(['mkdir', '-p', output_folder])

print 'counting how often every miRNA is gained or lost in every gene, in simulated and real datasets'
subprocess.call(['python', 'gene_specific_counts.py', simulation_folder, output_folder, miRNA_file])

print 'calculating simulation means and difference from simulated values for each gene in each miRNA and shuffle, and making output summary statistics'
subprocess.call(['python', 'analyze_single_removals_real_vs_fake2.py', output_folder, miRNA_file, master_folder])

print 'graphing'
subprocess.call(['python', 'graph_real_vs_fake.py', output_folder])
