'''
Prints all conserved miRNA seeds, their std dev values relative to simulated
values, and the miRNA names associated.

input files:
list of 8mers, mirnames (ex. real_eightmers_human_mouse_fish_better)
output files: real_beta_miRNA_(gains|losses|sum)
input scripts: populate_real_mers.py(AMTA stuff), make_distributions.py
output scripts: bin_distributions.py
'''
import custom
import sys
from decimal import Decimal
stat_folder=sys.argv[1]
stat_list=['gain', 'loss', 'sum']
#data_folder='/media/alfred/data_drive/big_data/turnover_rates/'
#conserved_file=open(data_folder+'real_eightmers_human_mouse_fish_better')
#conserved_file=open(data_folder+'real_eightmers_human_mouse_cow_not_fish')
conserved_file=open(sys.argv[2])

miRNA_list=[line.strip().split('\t') for line in conserved_file]
#print miRNA_list
conserved_miRNAs={line[0]:'['+line[1]+']' for line in miRNA_list}

for stat in stat_list:
	output2=open(stat_folder+'real_miRNA_'+stat, 'w')
	for line in open(stat_folder+'real_'+stat):
		line=line.strip().split('\t')
		eightmer, std_devs=line[0], line[1]
		if eightmer in conserved_miRNAs:
			output2.write('\t'.join([eightmer, std_devs, conserved_miRNAs[eightmer]])+'\n')
	output2.close()
