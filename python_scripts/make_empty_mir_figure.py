'''
Totally different from make_final_gene_lists2.py

Takes output from make_x_y.py and get_mir_exp.py. For each miRNA, returns the
tissue of highest expression, and the x and y coordinates for each gene in that
tissue (where x is absolute and y is relative expression of the gene).
'''
import cPickle
import random
import sys
data_folder=sys.argv[1]
empty_table={}

shared_tissues=dict([line.strip().split('\t') for line in open(data_folder+'tissue_key')][1:])

x_y_dict=cPickle.load(open(data_folder+'human_gene_atlas/x_and_y_dict', 'rb'))
mir_dict=cPickle.load(open(data_folder+'miRNA_expression/mir_gene_exp_dict3', 'rb'))

def remove_extra_tissues(input_dict, good_tissues):
	for gene in input_dict:
		for tissue in input_dict[gene][:]:
			if tissue[1] not in good_tissues:
				input_dict[gene].remove(tissue)
	return input_dict

shared_mir_dict=remove_extra_tissues(mir_dict, set(shared_tissues.keys()))

for mir in shared_mir_dict:
#	print mir_dict[mir]
	highest_value, highest_tissue=mir_dict[mir][-1]
	if highest_value>0:
		empty_table[mir]=[highest_tissue, x_y_dict[shared_tissues[highest_tissue]]]
cPickle.dump(empty_table, open(data_folder+'empty_figures', 'wb'), -1)
print 'done making empty mir figure'
