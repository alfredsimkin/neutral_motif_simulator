import cPickle
import custom
import copy
import time
import sys

data_folder, good_species=sys.argv[1], eval(sys.argv[2])
ancestor_dict=cPickle.load(open(data_folder+'ancestor_dict'))
gene_dict=cPickle.load(open(data_folder+'anc_aln'))
start_time=time.time()
new_dict={}

for gene in gene_dict:
	new_dict[gene]={}
	ref_seq=gene_dict[gene][good_species[0]]
	for bp_number, bp in enumerate(ref_seq):
		if bp=='N':
			for ancestor in ancestor_dict.values():
				unmodified_string=gene_dict[gene][ancestor]
				gene_dict[gene][ancestor]=unmodified_string[:bp_number]+'N'+unmodified_string[bp_number+1:]
	for descendant in ancestor_dict:
		ancestor=ancestor_dict[descendant]
		if ancestor not in new_dict[gene]:
			new_dict[gene][ancestor]={}
		new_dict[gene][ancestor][descendant]=[gene_dict[gene][ancestor].upper(), gene_dict[gene][descendant].upper()]
cPickle.dump(new_dict, open(data_folder+'reformatted_anc_dictU', 'wb'), -1)
cPickle.dump(new_dict.keys(), open(data_folder+'all_gene_names', 'wb'))
