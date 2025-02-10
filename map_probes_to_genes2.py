'''
Gathers statistics on affymetrix annotations and refseq genes, maps these to my
original list of refseq genes from my 3' UTR simulations, and removes probes
with _x_at in them (multiple mapping according to Farh et al. 2005)

This version(2) outputs a dictionary of genes with their associated probes, a
dictionary of probes with associated genes, and a dictionary of expression
levels for each probe
'''
import cPickle
import sys
data_folder, affy_ann_file, gnf1h_ann_file, exp_file, simulation_file=sys.argv[1:]
def map_probes(annotation_file, probe_column, gene_column, probe_dict, gene_dict, sim_genes):
	'''
	make a probe_dict with probe keys and refseq values, and a gene_dict with
	refseq keys and probe values
	'''
	no_title=True
	for line in open(annotation_file):
		if line.startswith('#'):
			continue
		if no_title:
			title=line.strip().split('\t')
			no_title=False
		else:
			split_line=line.strip('\n').split('\t')
			probe, gene_list=split_line[probe_column], split_line[gene_column]
			if '_x_at' not in probe:
				gene_list=gene_list.split(' /// ')
				for gene in gene_list:
					if gene in sim_genes:
						if gene not in gene_dict:
							gene_dict[gene]=[]
						gene_dict[gene].append(probe)
						if probe not in probe_dict:
							probe_dict[probe]=[]
						probe_dict[probe].append(gene)
	return probe_dict, gene_dict


def make_count_dict(input_dict):
	count_dict={}
	for the_key in input_dict:
		count=len(input_dict[the_key])
		if count not in count_dict:
			count_dict[count]=0
		count_dict[count]+=1
	return count_dict

def get_probe_exp(exp_file):
	'''
	takes each probe and puts tissue_specific expression levels in a dict
	'''
	exp_dict={}
	exp_table=[line.strip().split(',') for line in open(exp_file)]
	tissues=exp_table[0]
	exp_table=exp_table[1:]
	for line in exp_table:
		probe=line[0]
		if probe not in exp_dict:
			exp_dict[probe]={}
			for tissue_number, tissue in enumerate(tissues[1:]):
				if tissue not in exp_dict[probe]:
					exp_dict[probe][tissue]=float(line[tissue_number+1])
	return exp_dict

sim_genes=set(cPickle.load(open(simulation_file, 'rb')).keys())
exp_dict=get_probe_exp(exp_file)
probe_dict, gene_dict={},{}
probe_dict, gene_dict=map_probes(affy_ann_file, 0, 23, probe_dict, gene_dict, sim_genes)
probe_dict, gene_dict=map_probes(gnf1h_ann_file, 0, 2, probe_dict, gene_dict, sim_genes)
probe_counts=make_count_dict(probe_dict)
gene_counts=make_count_dict(gene_dict)
print 'probe counts:', probe_counts
print 'gene counts:', gene_counts
annot_gene_set=set(gene_dict.keys())
annot_probe_set=set(probe_dict.keys())

cPickle.dump(exp_dict, open(data_folder+'exp_dict', 'wb'), -1)
cPickle.dump(probe_dict, open(data_folder+'probe_dict', 'wb'), -1)
cPickle.dump(gene_dict, open(data_folder+'gene_dict', 'wb'), -1)
