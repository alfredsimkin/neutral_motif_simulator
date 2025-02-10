'''
calculates the proportion of miRNAs lying at a given number of std devs, as
well as the proportion of all 8mers lying at this point, and the binomial
probability of observing the proportion of miRNAs given the underlying
proportion of all 8mers.

input files: real_(gain, loss, sum), mir_ranks_(gain, loss, sum)
output files: distribution_(gain, loss, sum)
input scripts: make_distributions.py, extract_miRNA_ranks.py
'''
def make_list_dict(step_size, minimum, maximum):
	list_thing, dict_thing=[],{}
	counter=Decimal(minimum)
	while counter<=maximum:
		list_thing.append(counter)
		dict_thing[counter]=0
		counter+=step_size
	return list_thing, dict_thing

def map_list(value, step_size, list_thing, dict_thing):
	add_factor=step_size/Decimal(2)
	for list_item in list_thing:
		if value>(list_item-add_factor) and value<(list_item+add_factor):
			dict_thing[list_item]+=1
	return dict_thing


def print_probabilities(real_dict, full_dict, total_miRNAs, total_mers, x_range):
	output_file.write('rank\treal_successes\tfull_successes\treal_prob\tfull_prob\tunexpected\tcolor\n')
	for key in x_range:
		if key in real_dict:
			real_success=real_dict[key]
		else:
			real_success=0
		if key in full_dict:
			full_success=full_dict[key]
		else:
			full_success=0
		real_prob=real_success/float(total_miRNAs)
		full_prob=full_success/float(total_mers)
		cum_prob=Decimal(0)
		if real_prob<full_prob:
			unexpected=Decimal(1)-custom.cum_binom_prob(full_success, total_mers, real_success+1, real_count)
		else:
			unexpected=custom.cum_binom_prob(full_success, total_mers, real_success, real_count)
		if unexpected<0.001:
			color='red'
		elif unexpected<0.01:
			color='orange'
		elif unexpected<0.05:
			color='DarkViolet'
		else:
			color='black'
		output_file.write('\t'.join(map(str, [key, real_success, full_success, real_prob, full_prob, unexpected, color]))+'\n')

from decimal import Decimal
import custom
import math
import subprocess
import sys
steps=60
start=-30
end=30
step=Decimal(end-start)/Decimal(steps)
graph_types=['gain', 'loss', 'sum']
stats_folder=sys.argv[1]

for graph_type in graph_types:
	print graph_type
	output_file=open(stats_folder+'distribution_'+graph_type, 'w')
	real_file=open(stats_folder+'mir_ranks_'+graph_type)
	full_file=open(stats_folder+'real_'+graph_type)
	real_sum_dict, full_sum_dict={},{}
	full_count, real_count=0,0
	junk, real_sum_dict=make_list_dict(step, start, end)
	sig_list, full_sum_dict=make_list_dict(step, start, end)
	for line in real_file:
		real_count+=1
		sig_value=Decimal(line.strip().split('\t')[2])
		real_sum_dict=map_list(sig_value, step, sig_list, real_sum_dict)
	for line_number, line in enumerate(full_file):
		if line_number%1000==0:
			print line_number/65536.0
		full_count+=1
		sig_value=Decimal(line.strip().split('\t')[1])
		full_sum_dict=map_list(sig_value, step, sig_list, full_sum_dict)
	print_probabilities(real_sum_dict, full_sum_dict, real_count, full_count, sig_list)
