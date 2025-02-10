'''
unlike make_distributions.py, this script only rearranges my dataset of 1,000
dictionaries, each containing 65,536 8mers, each containing 3 values, into one
dictionary of 8mers, each 8mer containing 3 values (sites/gains/losses) each
containing 1,000 values (replicates).
'''
import custom
import cPickle
import subprocess
import sys
full_dict={}
data_folder='/media/alfred/data_drive/big_data/turnover_rates/'
sim_type=sys.argv[1]
replicate_start=int(sys.argv[2])
replicate_end=int(sys.argv[3])
folder=data_folder+sim_type+'/simulations/'
for replicate in range(replicate_start, replicate_end):
	print replicate
	input_dict=cPickle.load(open(folder+'eightmer_counts'+str(replicate)))
	for eightmer in input_dict:
		if 'N' in eightmer:
			print eightmer, 'is weird!!! stop!'
		if eightmer not in full_dict:
			full_dict[eightmer]={3:[]}
		for column_number, column in enumerate(input_dict[eightmer]):
			if column_number not in full_dict[eightmer]:
				full_dict[eightmer][column_number]=[]
			full_dict[eightmer][column_number].append(input_dict[eightmer][column_number])
		full_dict[eightmer][3].append(input_dict[eightmer][1]+input_dict[eightmer][2])
cPickle.dump(full_dict, open(data_folder+sim_type+'/full_raw_counts.cPickle', 'wb'), -1)
