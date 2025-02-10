'''
locates 8mers in high (or low) AU contexts in an alignment
'''
import cPickle
import sys

data_folder, high_low, mask_file=sys.argv[1:]
gene_dict=cPickle.load(open(data_folder+'reformatted_anc_dictU', 'rb'))
def modify_counts(anc_seq, desc_seq, high_low):
	sites=set([])
	for bp_number, anc_bp in enumerate(anc_seq):
		AU=0
		before, after=anc_seq[bp_number-30:bp_number], anc_seq[bp_number+8:bp_number+38]
		AU=before.count('A')+before.count('T')
		AU+=after.count('A')+after.count('T')
		if high_low=='low':
			if AU<(len(before)+len(after))/2.0:
				sites.add(bp_number)
		elif high_low=='high':
			if AU>(len(before)+len(after))/2.0:
				sites.add(bp_number)
	return sites

output_sites={}
for gene_number, gene in enumerate(gene_dict):
	output_sites[gene]={}
	for ancestor in gene_dict[gene]:
		output_sites[gene][ancestor]={}
		for descendant in gene_dict[gene][ancestor]:
			anc_seq, desc_seq=gene_dict[gene][ancestor][descendant]
			sites=modify_counts(anc_seq, desc_seq, high_low)
			output_sites[gene][ancestor][descendant]=sites
	if gene_number%100==0:
		print gene_number/float(len(gene_dict)), len(anc_seq), len(output_sites[gene]['6']['hg38']) #you will need to change this to match your reformatted_anc_dictU
cPickle.dump(output_sites, open(mask_file, 'wb'), -1)
