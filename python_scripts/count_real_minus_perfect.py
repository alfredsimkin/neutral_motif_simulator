'''
Counts the number of times a given 8mer occurs, is gained, and is lost in the
real dataset of 3' UTRs. Used to produce the 'observed' portion of observed vs.
expected turnover events.
'''
import custom
import cPickle
import sys
data_folder, stats_folder, mask_file=sys.argv[1:]
print 'loading mask'
mask_dict=cPickle.load(open(mask_file, 'rb'))
print 'mask loaded'

def modify_counts(eightmer_counts, anc_seq, desc_seq, mask_set):
	anc_seq=anc_seq.upper()
	desc_seq=desc_seq.upper()
	for bp_number, bp in enumerate(anc_seq):
		anc_8mer=anc_seq[bp_number:bp_number+8]
		eightmer_set=set(range(bp_number, bp_number+8))
		if len(eightmer_set&mask_set)==0:
			continue
#		if bp_number not in mask_set:
#			print anc_seq[bp_number-30:bp_number+38]
#			continue
		if len(anc_8mer)==8 and anc_8mer in good_set:
			desc_8mer=desc_seq[bp_number:bp_number+8]
			if desc_8mer in good_set:
				if anc_8mer not in eightmer_counts:
					eightmer_counts[anc_8mer]=[0,0,0]
				if desc_8mer not in eightmer_counts:
					eightmer_counts[desc_8mer]=[0,0,0]
				eightmer_counts[desc_8mer][0]+=1
				if desc_8mer!=anc_8mer:
					eightmer_counts[anc_8mer][2]+=1
					eightmer_counts[desc_8mer][1]+=1
	return eightmer_counts
eightmer_list=custom.count_in_base('AAAAAAAA', 4, 'ACGTz')
good_set=set(eightmer_list)
eightmer_counts={}
#gene_dict=cPickle.load(open(data_folder+'best_isoform_dictU', 'rb'))
gene_dict=cPickle.load(open(data_folder+'reformatted_anc_dictU', 'rb'))
sorted_genes=sorted(gene_dict.keys())
for gene in sorted_genes:
	for ancestor in gene_dict[gene]:
		for descendant in gene_dict[gene][ancestor]:
			anc, desc=gene_dict[gene][ancestor][descendant]
			eightmer_counts=modify_counts(eightmer_counts, anc, desc, mask_dict[gene])
cPickle.dump(eightmer_counts, open(stats_folder+'real_counts', 'wb'), -1)
#cPickle.dump(eightmer_counts, open(data_folder+'real_counts_best_isoform', 'wb'), -1)
