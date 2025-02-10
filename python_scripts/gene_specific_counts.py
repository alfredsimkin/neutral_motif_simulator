'''
A program that prints how many gains or losses are associated with every gene
in an alignment, for a small set of target eightmers
'''
import custom
import cPickle
import sys
input_data_folder, output_data_folder, mirfile=sys.argv[1:]
target_eightmers=set([line.strip().split('\t')[0] for line in open(mirfile)])

def modify_counts(eightmer_counts, anc_seq, desc_seq):
	anc_seq=anc_seq.upper()
	desc_seq=desc_seq.upper()
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
				if desc_8mer!=anc_8mer:
					eightmer_counts[anc_8mer][2]+=1
					eightmer_counts[desc_8mer][1]+=1
	return eightmer_counts

def count_alternate2(eightmer_counts, alternate_counts, gene):
	for eightmer in eightmer_counts:
		gain_count=eightmer_counts[eightmer][1]
		loss_count=eightmer_counts[eightmer][2]
		if (eightmer in target_eightmers) and (gain_count>0 or loss_count>0):
			if gene not in alternate_counts:
				alternate_counts[gene]={}
			alternate_counts[gene][eightmer]=[gain_count, loss_count]
			if eightmer not in alternate_counts['all']:
				alternate_counts['all'][eightmer]=[0,0]
			alternate_counts['all'][eightmer][0]+=gain_count
			alternate_counts['all'][eightmer][1]+=loss_count
	return alternate_counts

eightmer_list=custom.count_in_base('AAAAAAAA', 4, 'ACGTz')
good_set=set(eightmer_list)

for replicate in range(100):
	gene_dict=cPickle.load(open(input_data_folder+'/formatted_sim_wrong_anc_dict'+str(replicate), 'rb'))
	print 'replicate number is', replicate, 'out of 100 (hard coded)'
	print len(gene_dict)
	sorted_genes=sorted(gene_dict.keys())
	alternate_counts={'all':{}}
	for gene_number, gene in enumerate(sorted_genes):
		if gene_number%100==0:
			print gene_number
		eightmer_counts={}
		for ancestor in gene_dict[gene]:
			for descendant in gene_dict[gene][ancestor]:
				anc, desc=gene_dict[gene][ancestor][descendant]
				eightmer_counts=modify_counts(eightmer_counts, anc, desc)
		alternate_counts=count_alternate2(eightmer_counts, alternate_counts, gene)
	cPickle.dump(alternate_counts, open(output_data_folder+'miRNA_sim_counts'+str(replicate), 'w'))

gene_dict=cPickle.load(open(input_data_folder+'../reformatted_anc_dictU', 'rb'))
print len(gene_dict)
sorted_genes=sorted(gene_dict.keys())
alternate_counts={'all':{}}
for gene_number, gene in enumerate(sorted_genes):
	if gene_number%100==0:
		print gene_number
	eightmer_counts={}
	for ancestor in gene_dict[gene]:
		for descendant in gene_dict[gene][ancestor]:
			anc, desc=gene_dict[gene][ancestor][descendant]
			eightmer_counts=modify_counts(eightmer_counts, anc, desc)
	alternate_counts=count_alternate2(eightmer_counts, alternate_counts, gene)
cPickle.dump(alternate_counts, open(output_data_folder+'miRNA_gene_counts', 'w'))
'''
some scrap code, kept just in case:
def count_alternate(eightmer_counts, alternate_counts, gene):
	for eightmer in eightmer_counts:
		gain_count=eightmer_counts[eightmer][1]
		loss_count=eightmer_counts[eightmer][2]
		if (eightmer in target_eightmers) and (gain_count>0 or loss_count>0):
			if eightmer not in alternate_counts:
				alternate_counts[eightmer]={'gain':{}, 'loss':{}}
			if eightmer_counts[eightmer][1] not in alternate_counts[eightmer]['gain']:
				alternate_counts[eightmer]['gain'][gain_count]=set([])
			if eightmer_counts[eightmer][2] not in alternate_counts[eightmer]['loss']:
				alternate_counts[eightmer]['loss'][loss_count]=set([])
			alternate_counts[eightmer]['gain'][gain_count].add(gene)
			alternate_counts[eightmer]['loss'][loss_count].add(gene)
	return alternate_counts
'''
