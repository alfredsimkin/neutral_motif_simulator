'''
Used to make Table S15. This is a radically simplified version of
get_mir_exp3.py that uses miRNA seeds instead of attempting to standardize
miRNA names. It produces results identical to those when get_mir_exp3.py is
used (but with a radically different output file format geared specifically to
table S15).
'''
import cPickle
import sys
import custom

root_folder='/home/alfred/Dropbox/Simkin_Grimson_final_ms_folder/figures_and_data/mammal_data/deeply_conserved_mirs(human_mouse_fish)/'
exp_folder=root_folder+'coexpression_data/miRNA_expression/'
loss_file=root_folder+'main_results(vanilla)/mir_ranks_loss'
gain_file=root_folder+'main_results(vanilla)/mir_ranks_gain'
output_file=open(exp_folder+'table_S15.tsv', 'w')

mir_list=[line.strip().split('\t') for line in open(loss_file)] #table S1 of this work
input_exp_table=[line.strip().split('\t') for line in open(exp_folder+'human_w_random_precursor.tsv')] #Table S9 downloaded from Landgraf et al. 2007

def make_mir_exp_dict(input_exp_table):
	'''
	returns a dictionary of the fraction of all miRNA clones found in a
	tissue deriving from any particular miRNA seed (reverse complemented as an
	8mer target site
	'''
	exp_dict={}
	tissues=input_exp_table[0][8:] #names of all tissues (row 0 is title row)
	total_dict={}
	input_exp_table=input_exp_table[1:] #remove title row
	for line in input_exp_table:
		for tissue_number, tissue in enumerate(tissues):
			if tissue not in total_dict:
				total_dict[tissue]=0
			total_dict[tissue]+=float(line[tissue_number+8]) #how many miRNA clones total were recovered for each tissue?
	for line in input_exp_table:
		mirseed=custom.revcom(line[4][1:8])+'A'
		if mirseed not in exp_dict:
			exp_dict[mirseed]={}
		for tissue_number, tissue in enumerate(tissues):
			if tissue not in exp_dict[mirseed]:
				exp_dict[mirseed][tissue]=0
			exp_dict[mirseed][tissue]+=float(line[tissue_number+8])/total_dict[tissue] #what fraction of all the miRNA clones from this tissue derive from this seed?
	return exp_dict

def get_group(gain_rank, loss_rank):
	if gain_rank<3277 and loss_rank<3277:
		group=1
	elif gain_rank>=3277 and loss_rank<3277:
		group=2
	elif gain_rank>=3277 and loss_rank>=3277:
		group=3
	elif gain_rank<3277 and loss_rank>=3277:
		group=4
	return group

def get_exp_stats(tissue_dict):
	'''
	gets expression statistics from the 
	'''
	any_exp,above_1,max_exp=0,0,0
	tissue_list=[tissue_dict[tissue] for tissue in tissue_dict]
	for fraction in tissue_list:
		if fraction>0:
			any_exp+=1
		if fraction>0.01:
			above_1+=1
		if fraction>max_exp:
			max_exp=fraction
	median_exp=custom.get_median(tissue_list)
	avg_exp, var=custom.meanvar(tissue_list)
	return [any_exp, above_1, avg_exp, median_exp, max_exp]

def print_output_table(mir_exp_dict, loss_file, gain_file, output_file):
	'''
	extracts ranks from the gain and loss files and tissue expression
	statistics from mir_exp_dict and outputs these to a file
	'''
	gain_list=[line.strip().split('\t') for line in open(gain_file)]
	loss_list=[line.strip().split('\t') for line in open(loss_file)]
	mirseed, rank, name=0,1,5
	gain_dict={line[mirseed]:[line[rank], line[name]] for line in gain_list}
	loss_dict={line[mirseed]:[line[rank], line[name]] for line in loss_list}
	output_file.write('mirseed\tmir_name\tgain_rank\tloss_rank\tgroup\tany_exp\tabove_01\tAvg. % made up by this miRNA in any given tissue\tmedian\tmax_exp\n')
	for mirseed in gain_dict:
		mir_name=gain_dict[mirseed][1]
		gain_rank=int(gain_dict[mirseed][0])
		loss_rank=int(loss_dict[mirseed][0])
		group=get_group(gain_rank, loss_rank)
		output_list=map(str, [mirseed, mir_name, gain_rank, loss_rank, group])
		output_list+=map(str, get_exp_stats(mir_exp_dict[mirseed]))
		output_file.write('\t'.join(output_list)+'\n')

mir_exp_dict=make_mir_exp_dict(input_exp_table) #nested dictionary [broad_miRNA][narrow_miRNA][tissue] of fraction of all the miRNA clones from a given tissue that derive from a narrow miRNA copy
print mir_exp_dict['ACCAAAGA']
print_output_table(mir_exp_dict, loss_file, gain_file, output_file)
