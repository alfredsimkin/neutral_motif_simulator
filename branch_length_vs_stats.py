'''
plots branch length against some stats for the deeply conserved miRNA subset,
including overall deficit in number of binding sites and overall rank sum.
'''
import custom
import ete2
from ete2 import Tree
import cPickle
from decimal import Decimal
#phylogeny='/mnt/hulk_tank/alfred_stuff/silver_surfer/Grimson_new_phylogeny/scripts/data_prep_scripts_paper/outtree_labeled.nh'
phylogeny='/mnt/hulk_tank/alfred_stuff/silver_surfer/Grimson_new_phylogeny/data/other_data/branch_specific_stats/intree3.nh'
gene_dict=cPickle.load(open('/mnt/hulk_tank/alfred_stuff/silver_surfer/Grimson_new_phylogeny/data/other_data/reformatted_anc_dictU', 'rb'))
branch_folder='/mnt/hulk_tank/alfred_stuff/silver_surfer/Grimson_new_phylogeny/data/other_data/branch_specific_stats/'
gene=gene_dict.keys()[0]

pairing_list=[]
for ancestor in gene_dict[gene]:
	for descendant in gene_dict[gene][ancestor]:
		pairing_list.append((ancestor, descendant))

#pairing_list.remove(('1', 'monDom5'))
#pairing_list.remove(('3', 'mm10'))
t=Tree(phylogeny, format=1)
length_dict={}
for pair in pairing_list:
	descendant=pair[1]
	print descendant
	branch=t.search_nodes(name=descendant)[0]
	branch_length=Decimal(branch.dist)
	print branch_length
	length_dict[pair]=branch_length

for stat in ['gain', 'loss']:
	rank_sums=[]
	lengths=[]
	deficits=[]
	std_dev_sums=[]
	for pair in pairing_list:
		ancestor, descendant=pair
		stringed=ancestor+'-'+descendant
		lengths.append(float(length_dict[pair]))
		rank_sums.append(sum([int(line.split('\t')[1]) for line in open(branch_folder+stringed+'_mir_ranks_'+stat)]))
		observed_sum=(sum([float(line.split('\t')[3]) for line in open(branch_folder+stringed+'_mir_ranks_'+stat)]))
		expected_sum=(sum([float(line.split('\t')[4]) for line in open(branch_folder+stringed+'_mir_ranks_'+stat)]))
		std_dev_sums.append(sum([float(line.split('\t')[2]) for line in open(branch_folder+stringed+'_mir_ranks_'+stat)]))
		print expected_sum-observed_sum, pair, length_dict[pair]
		deficits.append(float(expected_sum-observed_sum))
	print lengths, rank_sums, deficits
	custom.scatter_it(lengths, rank_sums, 'branch_length', 'rank_sums_'+stat, faded=False)
	custom.scatter_it(lengths, deficits, 'branch_length', 'site_deficits_'+stat, lims=[[0,0.35], [0,2200]], faded=False)
#	custom.scatter_it(lengths, deficits, 'branch_length', 'site_deficits_'+stat, lims=[[0,0.35], [-400,820]], faded=False)
	custom.scatter_it(lengths, std_dev_sums, 'branch_length', 'std_dev_sums_'+stat, faded=False)
