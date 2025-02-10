'''
goes through a list of all species in mirbase 21 and crossreferences a 100way
phylogeny from UCSC hg19 (29 sp overlap). Lists which species from the
overlapping list a full-length mature miRNA is found in, and calculates a
branch sum over which each full-length miRNA is conserved. Outputs each full
mature miRNA, 8mer, and branch length sum.

Dependencies:
1. ete2 python package in a unix $PYTHONPATH folder (from etetoolkit.org)
2. tree_manips.py (provided here)
3. custom.py (provided here)
4. mirbase mature sequences (downloaded from mirbase)
5. mirbase full scientific names (downloaded from mirbase)
5. phylogeny relating mirbase species to each other (downloaded from UCSC)

Input parameters:
1. name of mirbase mature sequence file
2. name of mirbase species scientific names file
3. name of newick formatted phylogeny file
4. desired name of output newick phylogeny (which names internal nodes)
5. name of desired output file with conservation levels of each mature miRNA
'''
import tree_manips
import custom

#input parameters below:
mature_file='mirbase_21_mature'
mirbase='organisms.txt'
hundred_way='hg38.100way.scientificNames.nh'
outtree=hundred_way[:-3]+'mod.nh'
output_file=open('mature_mir_branch_lengths', 'w')

def get_overlap(mirbase, hundred_way):
	mirbase=[line.strip().split('\t') for line in open(mirbase)]
	species_dict={}
	for line in mirbase:
		species=line[2].replace(' ', '_')
		short=line[0]
		for line in open(hundred_way):
			if species in line:
				species_dict[short]=species
	return species_dict

def name_internal(input_file):
	string, new_string='',''
	for line in open(input_file):
		string+=line.strip()
	string_list=string.split('):')
	for counter, item in enumerate(string_list[:-1]):
		new_string+=item+')'+str(counter)+':'
	new_string+=string_list[-1]
	return new_string

def parse_mature(mirbase_table, overlapping_species):
	new_dict={}
	for line in open(mirbase_table):
		line=line.strip().split('\t')
		short=line[0][:3]
		if short in overlapping_species:
			if line[1] not in new_dict:
				new_dict[line[1]]=[]
			if overlapping_species[short] not in new_dict[line[1]]:
				new_dict[line[1]].append(overlapping_species[short])
	return new_dict

overlapping_species=get_overlap(mirbase, hundred_way)
print overlapping_species
named_tree=name_internal(hundred_way)
output_newick=open(outtree, 'w')
output_newick.write(named_tree)
output_newick.close()
mature_dict=parse_mature(mature_file, overlapping_species)
for mature in mature_dict:
	if len(mature_dict[mature])>1:
		stringed=map(str, ([mature, custom.revcom(mature[1:8])+'A', tree_manips.get_dist(outtree, mature_dict[mature])]))
		output_file.write('\t'.join(stringed)+'\n')
