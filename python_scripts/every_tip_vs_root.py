'''
Takes existing simulations and makes a new ancestor_descendant alignment that
compares every tip species to the common root
'''
root_node='1'
species_list=['hg19', 'nomLeu3', 'calJac3', 'otoGar3', 'canFam3', 'felCat5', 'equCab2', 'oviAri3', 'bosTau7', 'mm10']
real_path='/home/alfred/grimson_postdoc/big_data/turnover_rates/4-28-15_deeper_phylogeny/reformatted_anc_dictU'
sim_path='/home/alfred/grimson_postdoc/big_data/turnover_rates/4-28-15_deeper_phylogeny/simulations/'

def revise_dict(input_dict):
	for gene in input_dict:
		root_seq=input_dict[gene][root_node][0]
		if gene not in new_dict:
			new_dict[gene]={root:{}}
		for ancestor in input_dict[gene]:
			for descendant in input_dict[gene][ancestor]:
				descendant_seq=input_dict[gene][ancestor][descendant][1]
				if descendant in species_list:
					new_dict[root_node][descendant]=[root_seq, descendant_seq]
	return new_dict

real_input=cPickle.load(open(real_path+'reformatted_anc_dictU', 'rb'))
revised_real=revise_dict(real_input)
cPickle.dump(revised_real, open(real_path+'reformatted_anc_dictU_root', 'wb'))

for rep in range(100):
	sim_input=cPickle.load(open(sim_path+'formatted_sim_wrong_anc_dict'+str(rep), 'rb'))
	revised_sim=revise_dict(sim_input)
	cPickle.dump(revised_sim, open(sim_path+'formatted_sim_anc_dict_root'+str(rep), 'wb'))
