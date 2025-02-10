'''
Replaces get_mir_exp2.py

calculates human miRNA expression levels using clones mapping to the 5p or 3p
arm. Clones mapping to multiple precursors are noted. Often multiple miRNAs
share the same seed. In these cases, I'll take the sum. Unlike get_mir_exp.py,
returns mir expression in each tissue as a fraction of all clones in the tissue
(summed with total_dict) rather than as a raw number

This script could be radically simplified by reverse complementing the
seed to match the 'most matched sequence' column of 'human_w_random_
precursor.tsv' (an approach used by get_mir_exp4.py to produce table S15)

Called by get_coexpression.py
'''
import cPickle
import sys

data_folder, mir_file, full_mir_file=sys.argv[1:]
mir_list=[line.strip().split('\t') for line in open(mir_file)] #mir_file is 'mir_ranks_loss' from main pipeline
full_mir_table=[line.strip().split('\t') for line in open(full_mir_file)] #full_mir_file is mirbase_21_table_fixed.tsv
manual_mappings={'hsa-miR-451a':'5p', 'hsa-miR-217':'5p'}
#mir_exp_table=[line.strip().split('\t') for line in open(data_folder+'miRNA_expression/human_w_ambig_precursor.tsv')]
mir_exp_table=[line.strip().split('\t') for line in open(data_folder+'human_w_random_precursor.tsv')] #downloaded from Landgraf et al. 2007

def get_mir_info(mir):
	'''
	crudely maps mirs without arm info to the 5p or 3p arm by cutting the
	precursor approximately in half. Two mirs are in neither half and therefore
	manually mapped
	'''
	found=False
	if mir in manual_mappings:
		return manual_mappings[mir]
	else:
		for line in full_mir_table:
			if line[4]==mir:
				found=True
				cut_point=len(line[2])/2
				if line[5] in line[2][:cut_point]:
					result='5p'
				elif line[5] in line[2][cut_point:]:
					result='3p'
				else:
					print 'snapped'
		if not found:
			print 'Error: missing'
		return result

def make_mir_exp_dict(mir_exp_table):
	'''
	returns a nested dictionary of the fraction of all miRNA clones found in a
	tissue deriving from any particular miRNA
	'''
	exp_dict={}
	tissues=mir_exp_table[0][8:] #names of all tissues (row 0 is title row)
	total_dict={}
	mir_exp_table=mir_exp_table[1:] #remove title row
	for line in mir_exp_table:
		for tissue_number, tissue in enumerate(tissues):
			if tissue not in total_dict:
				total_dict[tissue]=0
			total_dict[tissue]+=float(line[tissue_number+8]) #how many miRNA clones total were recovered for each tissue?
	for line in mir_exp_table:
		mirs=line[0].split(', ') #this line probably isn't needed
		for mir in mirs:
			new_mir='-'.join(mir.split('-')[:3]) #remove identical mir copies from names (e.g. hsa-mir-92a-1 becomes hsa-mir-92a)
			if new_mir[-1].isalpha(): #remove nonidentical mir copies from names (e.g. hsa-mir-92a becomes hsa-mir-92)
				new_mir=new_mir[:-1]
			new_mir=new_mir+'-'+line[2] #add '5p' or '3p' to mir name
			if not line[0].startswith('hsa') or int(float(line[7]))==0: #ignore non-human miRNAs
				continue
			if new_mir not in exp_dict:
				exp_dict[new_mir]={}
			if mir not in exp_dict[new_mir]:
				exp_dict[new_mir][mir]={}
			for tissue_number, tissue in enumerate(tissues):
				if tissue not in exp_dict[new_mir][mir]:
					exp_dict[new_mir][mir][tissue]=float(line[tissue_number+8])/total_dict[tissue] #what fraction of all the miRNA clones from this tissue derive from this miRNA?
#	for sub_mir in exp_dict['hsa-mir-302-3p']:
#		print sub_mir, exp_dict['hsa-mir-302-3p'][sub_mir]['hsa_Seminoma-2073']#['hsa-mir-302c']['hsa_Seminoma-2073']
	return exp_dict

def get_tissue_values(mir_dict, tissue_values):
	'''
	returns a tissue-centric dictionary containing a list of the fractions of
	all clones for a given tissue derive from all copies of a particular miRNA
	(e.g. mir-92a-1+mir-92a-2+mir-92b). Called multiple times to retrieve
	values across multiple miRNAs that happen to share a seed (e.g.
	hsa-miR-92a-3p and hsa-miR-363-3p)
	'''
	for mir in mir_dict:
		for tissue in mir_dict[mir]:
			if tissue not in tissue_values:
				tissue_values[tissue]=[]
			tissue_values[tissue].append(float(mir_dict[mir][tissue]))
	return tissue_values
	
def get_sum(tissue_values):
	for tissue in tissue_values:
		tissue_values[tissue]=sum(tissue_values[tissue])
	return tissue_values

mir_exp_dict=make_mir_exp_dict(mir_exp_table) #nested dictionary [broad_miRNA][narrow_miRNA][tissue] of fraction of all the miRNA clones from a given tissue that derive from a narrow miRNA copy
summed_dict, mir_centric_mir, tissue_centric_mir={}, {}, {}
for mir_seed in mir_list:
	#('AGCACTTA', 'hsa-miR-302a-3p, hsa-miR-302b-3p, hsa-miR-302c-3p, hsa-miR-302d-3p, hsa-miR-372-3p, hsa-miR-373-3p, hsa-miR-520e, hsa-miR-520a-3p, hsa-miR-520b, hsa-miR-520c-3p, hsa-miR-520d-3p, hsa-miR-302e')
	mir_names=mir_seed[5].replace("'", '') #not needed
	new_mir_names=mir_names.split(', ')
#	print mir_names
	tissue_values={}
	seen_already=set([])
	for mir in new_mir_names:
		if mir.endswith('5p'):
			arm='5p'
		elif mir.endswith('3p'):
			arm='3p'
		else:
			arm=get_mir_info(mir) #tries to retrieve which arm of the precursor the mature miRNA derives from
		mir=mir.replace('-5p', '').replace('-3p', '').lower()
		if mir[-1].isalpha():
			mir=mir[:-1]
		mir=mir+'-'+arm #mir names from mir_ranks_loss (mirbase derived) should now be formatted identically to modified names from Landgraf et al
		if mir in mir_exp_dict and mir not in seen_already:
			tissue_values=get_tissue_values(mir_exp_dict[mir], tissue_values)
			seen_already.add(mir)
#		else:
#			tissue_values='not found'
	the_sum=get_sum(tissue_values)
	summed_dict[(mir_seed[0], mir_names)]=the_sum #returns a dictionary containing fraction of all miRNA clones in a given tissue that can be attributed to a given miRNA seed
#	if mir_seed[0]=='AGCACTTA':
#		print the_sum['hsa_Seminoma-2073']

for mir in summed_dict:
	if mir not in mir_centric_mir:
		mir_centric_mir[mir]=[]
#	if mir==('ACTTTATA', 'hsa-miR-142-5p'):
#	if mir==('GACAATCA', 'hsa-miR-219a-5p, hsa-miR-4782-3p, hsa-miR-6766-3p'):
#		print mir, summed_dict[mir]
	for tissue in summed_dict[mir]:
		if tissue not in tissue_centric_mir:
			tissue_centric_mir[tissue]=[]
		tissue_list=[summed_dict[mir][tissue], tissue]
		mir_list=[summed_dict[mir][tissue], mir]
		mir_centric_mir[mir].append(tissue_list)
		tissue_centric_mir[tissue].append(mir_list)
for tissue in tissue_centric_mir:
	tissue_centric_mir[tissue].sort()
for mir in mir_centric_mir:
	mir_centric_mir[mir].sort()
cPickle.dump(tissue_centric_mir, open(data_folder+'mir_tissue_exp_dict3', 'wb'), -1)
cPickle.dump(mir_centric_mir, open(data_folder+'mir_gene_exp_dict3', 'wb'), -1)
print 'done with get_mir_exp3.py'
#for mir in mir_centric_mir:
#	print mir, len(mir_centric_mir[mir])
