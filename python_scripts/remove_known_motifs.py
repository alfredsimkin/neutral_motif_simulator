'''
based off of make_control_kmers2.py

Written to look for putatively novel motifs. Removes any 8mers that correspond
to 6mer cores of known conserved miRNAs, PUF proteins (conservatively, ones
with TGT in them), TGCAT, AGGT, polyA sites (AATAAA) and low complexity motifs
(defined as those with 2 or fewer nucleotides).
'''
import custom
import cPickle
eightmer_list=custom.count_in_base('AAAAAAAA', 4, 'ACGTz')
mirseeds=[line.strip().split()[0] for line in open('merged_conserved_mirs')]
gain_list=[line.strip().split() for line in open('../figure_data_for_rene/mammalian_phylogeny_8mer_results/3UTR/real_gain')]
loss_list=[line.strip().split() for line in open('../figure_data_for_rene/mammalian_phylogeny_8mer_results/3UTR/real_loss')]
sub_size=6
subs=[[1,7]]
bad_kmers=['TGT', 'TGCAT', 'AGGT', 'AATAAA']
outfile=open('putative_novel_kmers', 'w')

def make_mirset(mirseeds, subs):
	'''
	lists the set of all 'core' subseeds of the mirseed (ex. nucleotides 1:7 as
	the core 6mer of the range 0:8)
	'''
	mirset=set([])
	for mirseed in mirseeds:
		for sub in subs:
			subseed=mirseed[sub[0]:sub[1]]
			if subseed not in mirset:
				mirset.add(subseed)
	return mirset

def remove_mirs(input_set, sub_size, mirset):
	'''
	looks at all windows of every 8mer and removes those matching a miRNA
	'core'
	'''
	import copy
	for eightmer in copy.deepcopy(input_set):
		for bp_number in xrange(8-sub_size+1):
			if eightmer[bp_number:bp_number+sub_size] in mirset:
				input_set.discard(eightmer)
	return list(input_set)

def remove_bad_kmers(good_list, bad_list):
	'''
	a fairly slow function that tests whether each 'bad' kmer is found anywhere
	within a larger 'good' kmer 
	'''
	for good in good_list[:]:
		for bad in bad_list:
			if bad in good and good in good_list:
				good_list.remove(good)
	return good_list

def remove_simple(good_list):
	'''
	removes any simple kmers from the list
	'''
	for good in good_list[:]:
		if len(set(good))<3:
			good_list.remove(good)
	return good_list

def print_sorted(good_eightmers, score_list, outfile_name):
	outfile=open(outfile_name, 'w')
	good_eightmers=set(good_eightmers)
	'''
	sorts all eightmers by their ranked stdev score, and prints those that are
	considered 'good' after the above filters.
	'''
	score_list=sorted([[float(line[1]), line[0]] for line in score_list])
	for line in score_list:
		if line[1] in good_eightmers:
			outfile.write(line[1]+'\t'+str(line[0])+'\n')

mirset=make_mirset(mirseeds, subs)
print 'removing mirs from eightmers'
eightmer_list=remove_mirs(set(eightmer_list), sub_size, mirset)
print 'removing other bad kmers from eightmers'
eightmer_list=remove_bad_kmers(eightmer_list, bad_kmers)
print 'removing simple kmers'
eightmer_list=remove_simple(eightmer_list)
print 'printing ranked gains'
print_sorted(eightmer_list, gain_list, 'good_gain_motifs')
print 'printing ranked losses'
print_sorted(eightmer_list, loss_list, 'good_loss_motifs')
