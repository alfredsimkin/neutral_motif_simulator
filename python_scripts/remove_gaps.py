'''
takes input dictionary from find_maf_overlaps.py, outputs a version with runs
of gaps and ambiguous nucleotides in one or more species replaced by a single
ambiguous nucleotide in all species
'''

import custom
import cPickle
import copy
import sys
data_folder, species_list, common_sp, prefix, input_dict=sys.argv[1:]
species_list, common_sp=eval(species_list), eval(common_sp)
strand_dict=dict([line.strip().split()[3:6:2] for line in open(data_folder+prefix+'.bed')])
gene_dict=cPickle.load(open(data_folder+input_dict, 'rb'))
good_set=set('ACGTacgt')
for gene_number, gene in enumerate(gene_dict):
	if gene_number%1000==0:
		print gene, float(gene_number)/len(gene_dict)
	unzipped_list=[]
	if len(gene_dict[gene])>1:
		print gene, 'is odd'
		print gene_dict[gene]
		exit()
	for coord in gene_dict[gene]:
		for species in species_list:
			unzipped_list.append(gene_dict[gene][coord][species])
		gene_dict[gene]={}
		zipped_list=[''.join(char) for char in zip(*unzipped_list)]
		new_zipped=[]
		for pos_number, pos in enumerate(zipped_list):
			if zipped_list[pos_number]=='-'*len(common_sp):
				pass
			elif len(set(zipped_list[pos_number])-good_set)>0:
#				print pos_number, len(zipped_list)
				if pos_number<len(zipped_list)-1 and len(set(zipped_list[pos_number+1])-good_set)==0:
					new_zipped.append('N'*len(species_list))
			else:
				new_zipped.append(zipped_list[pos_number])
		if new_zipped==[]:
			new_zipped.append('X'*len(species_list))
#		print zipped_list
#		print new_zipped
		if len(set(zipped_list[-1])-good_set)>0:
			new_zipped[-1]='N'*len(zipped_list[pos_number])
		new_alignment=[''.join(char) for char in zip(*new_zipped)]
#		print new_alignment
		for species_number, species in enumerate(species_list):
			common=common_sp[species_number]
			if strand_dict[gene]=='+':
				gene_dict[gene][common]=new_alignment[species_number].upper()
			elif strand_dict[gene]=='-':
				gene_dict[gene][common]=custom.revcom(new_alignment[species_number].upper())
			else:
				print 'ERROR!'
				exit()
for gene in gene_dict.keys():
	if len(gene_dict[gene][common_sp[0]])<50:
		gene_dict.pop(gene)
cPickle.dump(gene_dict, open(data_folder+input_dict+'_nogaps', 'wb'), -1)
