'''
Created to produce a more accurate model of dinucleotide turnover probabilities
given flanking nucleotide sequence, this program asks how often a particular
trinucleotide turns into another trinucleotide given that its first and third
nucleotides do not change. (i.e. how often a nucleotide changes given the
context it occurs in).
'''

import custom
import cPickle
import sys
kmer=1

#alignment_file='/media/alfred/data_drive/big_data/turnover_rates/1-27-15_nonredundant/simulations/formatted_sim_wrong_anc_dict0'
data_folder=sys.argv[1]
alignment_file=data_folder+'reformatted_anc_dictU'
def check_kmer(kmer_counts, anc_kmer, desc_kmer):
	'''
	configured for 1mers (monomers)
	'''
	if anc_kmer not in kmer_counts:
		kmer_counts[anc_kmer]={'change':0, anc_kmer:0}
	if desc_kmer not in kmer_counts[anc_kmer]:
		kmer_counts[anc_kmer][desc_kmer]=0
	if anc_kmer!=desc_kmer:
		kmer_counts[anc_kmer]['change']+=1
		kmer_counts[anc_kmer][desc_kmer]+=1
	else:
		kmer_counts[anc_kmer][anc_kmer]+=1
	return kmer_counts

def modify_counts(kmer_counts, anc_seq, desc_seq):
	anc_seq=anc_seq.upper()
	desc_seq=desc_seq.upper()
	for bp_number, bp in enumerate(anc_seq):
		anc_kmer=anc_seq[bp_number:bp_number+kmer]
		if len(anc_kmer)==kmer:
			desc_kmer=desc_seq[bp_number:bp_number+kmer]
			kmer_counts=check_kmer(kmer_counts, anc_kmer, desc_kmer)
	return kmer_counts
gene_dict=cPickle.load(open(alignment_file))
sorted_genes=sorted(gene_dict.keys())
kmer_counts={}
for gene in sorted_genes:
	for ancestor in gene_dict[gene]:
		if ancestor not in kmer_counts:
			kmer_counts[ancestor]={}
		for descendant in gene_dict[gene][ancestor]:
			if descendant not in kmer_counts[ancestor]:
				kmer_counts[ancestor][descendant]={}
			anc, desc=gene_dict[gene][ancestor][descendant]
			kmer_counts[ancestor][descendant]=modify_counts(kmer_counts[ancestor][descendant], anc, desc)
cPickle.dump(kmer_counts, open(data_folder+'real_'+str(kmer)+'mer_counts', 'wb'), -1)

'''
good_kmers=set(custom.count_in_base('AAA', 5, 'ACGTNz'))
maximum=1000
special_k=''
for kmer in kmer_counts['1']['2']:
	changers=kmer_counts['1']['2'][kmer]['change']
	not_changers=kmer_counts['1']['2'][kmer][kmer]
	print kmer, changers, not_changers
	if kmer in good_kmers and changers>0:
		current=not_changers/changers
		if current<maximum:
			maximum=current
			special_k=kmer
print '**', special_k, maximum
'''
