'''
Counts the number of times a given 8mer occurs, is gained, and is lost in the
real dataset of 3' UTRs. Used to produce the 'observed' portion of observed vs.
expected turnover events. Does this for every type of seed conversion event.
'''

def get_7mer1A(eightmer):
	eightmer_set=set([])
	constant=eightmer[1:]
	start_set=nucleotides-set(eightmer[0])
	for start in start_set:
		eightmer_set.add(start+constant)
	return eightmer_set
def get_6mer(eightmer):
	eightmer_set=set([])
	constant=eightmer[1:7]
	start_set=nucleotides-set(eightmer[0])
	end_set=nucleotides-set(eightmer[7])
	for start in start_set:
		for end in end_set:
			eightmer_set.add(start+constant+end)
	return eightmer_set
def get_7merm8(eightmer):
	eightmer_set=set([])
	constant=eightmer[:7]
	end_set=nucleotides-set(eightmer[7])
	for end in end_set:
		eightmer_set.add(constant+end)
	return eightmer_set

def make_events(eightmer_list):
	event_dict={}
	for eightmer in eightmer_list:
		event_dict[eightmer]={}
		event_dict[eightmer]['7mer1A']=get_7mer1A(eightmer)
		event_dict[eightmer]['7merm8']=get_7merm8(eightmer)
		event_dict[eightmer]['6mer']=get_6mer(eightmer)
		event_dict[eightmer]['8mer']=set([eightmer])
	return event_dict
def modify_counts(master_counts, anc_seq, desc_seq):
	anc_seq=anc_seq.upper()
	desc_seq=desc_seq.upper()
	for bp_number, bp in enumerate(anc_seq):
		anc_8mer=anc_seq[bp_number:bp_number+8]
		if len(anc_8mer)==8 and anc_8mer in good_set:
			desc_8mer=desc_seq[bp_number:bp_number+8]
			if desc_8mer in good_set:
				for seed_comp in master_counts:
					first_seed, second_seed=seed_comp.split('_')
					if anc_8mer not in master_counts[seed_comp]:
						master_counts[seed_comp][anc_8mer]=[0,0,0]
					if desc_8mer not in master_counts[seed_comp]:
						master_counts[seed_comp][desc_8mer]=[0,0,0]
					master_counts[seed_comp][anc_8mer][0]+=1
					master_counts[seed_comp][desc_8mer][0]+=1
					if desc_8mer!=anc_8mer:
						if first_seed!='none' and anc_8mer[1:7]==desc_8mer[1:7]:
							gain_8mers=event_dict[anc_8mer][first_seed]&event_dict[desc_8mer][second_seed]
							loss_8mers=event_dict[anc_8mer][second_seed]&event_dict[desc_8mer][first_seed]
						elif first_seed=='none' and anc_8mer[1:7]!=desc_8mer[1:7]:
							gain_8mers=event_dict[desc_8mer][second_seed]
							loss_8mers=event_dict[anc_8mer][second_seed]
						else:
							continue
						for gained_8mer in gain_8mers:
							if gained_8mer not in master_counts[seed_comp]:
								master_counts[seed_comp][gained_8mer]=[0,0,0]
							master_counts[seed_comp][gained_8mer][1]+=1
						for lost_8mer in loss_8mers:
							if lost_8mer not in master_counts[seed_comp]:
								master_counts[seed_comp][lost_8mer]=[0,0,0]
							master_counts[seed_comp][lost_8mer][2]+=1
	return master_counts

import custom
import cPickle
import sys
data_folder, stats_folder, seed_types=sys.argv[1:]
nucleotides=set('ACGT')
eightmer_list=custom.count_in_base('AAAAAAAA', 4, 'ACGTz')
good_set=set(eightmer_list)
seed_types=seed_types.split('_')

master_counts={}
for small_number, small_seed in enumerate(seed_types):
	for large_seed in seed_types[small_number+1:]:
		master_counts[small_seed+'_'+large_seed]={}

#gene_dict=cPickle.load(open(data_folder+'best_isoform_dictU', 'rb'))
gene_dict=cPickle.load(open(data_folder+'reformatted_anc_dictU', 'rb'))
sorted_genes=sorted(gene_dict.keys())
event_dict=make_events(eightmer_list)

for gene in sorted_genes:
	for ancestor in gene_dict[gene]:
		for descendant in gene_dict[gene][ancestor]:
			anc, desc=gene_dict[gene][ancestor][descendant]
			master_counts=modify_counts(master_counts, anc, desc)
cPickle.dump(master_counts, open(stats_folder+'everything_real_counts', 'wb'), -1)
#cPickle.dump(eightmer_counts, open(data_folder+'real_counts_best_isoform', 'wb'), -1)
