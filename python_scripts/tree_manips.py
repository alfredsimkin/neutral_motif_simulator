from ete2 import Tree
#t=Tree('100way_mod.nh', format=1)

def get_dist(t, species_list):
	from decimal import Decimal
	t=Tree(t, format=1)
	ancestor=t.get_common_ancestor(species_list)
	branch_set=set([])
	total_dist=Decimal('0')
	for species in species_list:
		branch=t.search_nodes(name=species)[0]
		while branch!=ancestor:
			branch_set.add(branch)
			branch=branch.up
	for branch in branch_set:
		total_dist+=Decimal(branch.dist)
	return total_dist

def prune_species(t, good_species, outname):
	t=Tree(t, format=1)
	t.prune(good_species, preserve_branch_length=True)
	t.write(format=1, outfile=outname)
