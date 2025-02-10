'''
This is the main script. It uses a table of substitution probabilities to
create random substitutions on a series of input ancestral sequences.
'''

import cPickle
import random
import custom
import copy
import subprocess
import bisect
import sys
data_folder=sys.argv[1]
gene_dict=cPickle.load(open(data_folder+'reformatted_anc_dictU', 'rb'))
sim_gene_dict=copy.deepcopy(gene_dict)
change_dict=cPickle.load(open(data_folder+'real_1mer_counts', 'rb'))
subprocess.call(['mkdir', data_folder+'simulations'])
good_set=set(custom.count_in_base('AAAAAAAA', 4, 'ACGTz'))
start_rep, end_rep=map(int, sys.argv[2:4])
reference, ref_anc=sys.argv[4:]

def create_lists(change_dict):
	list_dict=copy.deepcopy(change_dict)
	for ancestor in change_dict:
		for descendant in change_dict[ancestor]:
			for anc_1mer in change_dict[ancestor][descendant]:
				descendant_list=make_pair(change_dict[ancestor][descendant][anc_1mer])
				list_dict[ancestor][descendant][anc_1mer]=descendant_list
	return list_dict

def modify_counts(eightmer_counts, anc_seq, desc_seq):
	for bp_number, bp in enumerate(anc_seq):
		anc_8mer=anc_seq[bp_number:bp_number+8]
		if len(anc_8mer)==8 and anc_8mer in good_set:
			desc_8mer=desc_seq[bp_number:bp_number+8]
			if desc_8mer in good_set:
				if anc_8mer not in eightmer_counts:
					eightmer_counts[anc_8mer]=[0,0,0]
				if desc_8mer not in eightmer_counts:
					eightmer_counts[desc_8mer]=[0,0,0]
				eightmer_counts[desc_8mer][0]+=1
				if ancestor=='1':
					eightmer_counts[anc_8mer][0]+=1
				if desc_8mer!=anc_8mer:
					eightmer_counts[anc_8mer][2]+=1
					eightmer_counts[desc_8mer][1]+=1
	return eightmer_counts

def make_pair(input_dict):
	counter=0
	paired_list=[[],[]]
	for key in input_dict:
		if key!='change':
			paired_list[0].append(counter)
			paired_list[1].append(key)
			counter+=input_dict[key]
	paired_list[0].append(counter)
	return paired_list

def choice(paired_list):
	choice=random.randrange(paired_list[0][-1])
	index=bisect.bisect_right(paired_list[0], choice)
	seq=paired_list[1][index-1]
	return seq

def simulate_descendant(anc_seq, desc_seq, ancestor, descendant):
	sim_desc=anc_seq
	for bp_number, bp in enumerate(sim_desc[:-2]):
		anc_1mer=sim_desc[bp_number:bp_number+1]
		if anc_1mer in list_dict[ancestor][descendant]:
			descendant_list=list_dict[ancestor][descendant][anc_1mer]
			desc_1mer=choice(descendant_list)
			sim_desc=sim_desc[:bp_number]+desc_1mer+sim_desc[bp_number+1:]
	return sim_desc

list_dict=create_lists(change_dict)
for rep in range(start_rep, end_rep):
	eightmer_counts={}
	print 'rep is', rep
	for gene_number, gene in enumerate(gene_dict):
		if gene_number%100==0:
			print float(gene_number)/len(gene_dict)
		ref_seq=gene_dict[gene][ref_anc][reference][1]
		if len(ref_seq)<50:
			continue
		for ancestor in gene_dict[gene]:
			for descendant in gene_dict[gene][ancestor]:
				anc_seq=gene_dict[gene][ancestor][descendant][0]
				desc_seq=gene_dict[gene][ancestor][descendant][1]
				new_seq=simulate_descendant(anc_seq, desc_seq, ancestor, descendant)
				eightmer_counts=modify_counts(eightmer_counts, anc_seq, new_seq)
				sim_gene_dict[gene][ancestor][descendant][1]=new_seq
	cPickle.dump(sim_gene_dict, open(data_folder+'simulations/formatted_sim_wrong_anc_dict'+str(rep), 'wb'), -1)
	cPickle.dump(eightmer_counts, open(data_folder+'simulations/eightmer_counts'+str(rep), 'w'))
