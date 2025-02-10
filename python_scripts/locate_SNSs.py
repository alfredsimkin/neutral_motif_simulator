'''
locates single nucleotide substitutions in an input dataset
'''
import cPickle
import sys

data_folder, mask_file=sys.argv[1:]
gene_dict=cPickle.load(open(data_folder+'reformatted_anc_dictU', 'rb'))
def modify_counts(SNS_sites, anc_seq, desc_seq):
	for bp_number, anc_bp in enumerate(anc_seq):
		desc_bp=desc_seq[bp_number]
		if desc_bp!=anc_bp:
			SNS_sites.add(bp_number)
	return SNS_sites

output_sites={}
for gene in gene_dict:
	SNS_sites=set([])
	for ancestor in gene_dict[gene]:
		for descendant in gene_dict[gene][ancestor]:
			anc_seq, desc_seq=gene_dict[gene][ancestor][descendant]
			SNS_sites=modify_counts(SNS_sites, anc_seq, desc_seq)
	output_sites[gene]=SNS_sites
	print gene, len(anc_seq), len(output_sites[gene])
cPickle.dump(output_sites, open(mask_file, 'wb'), -1)
