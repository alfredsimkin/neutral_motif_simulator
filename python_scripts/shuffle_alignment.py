'''
shuffles alignments from 3' UTRs. Separates 3' UTRs into a first half and a
second half, and then interdigitates 5mer sites from the first half of the UTR
with 5mer sites from the second half of the UTR, such that no 5mer is adjacent
its original 5mer neighbor
'''
import copy
import cPickle

unshuffled_file='/mnt/hulk_tank/alfred_stuff/silver_surfer/Grimson_new_phylogeny/data/other_data/reformatted_anc_dictU'
shuffled_file='/mnt/hulk_tank/alfred_stuff/silver_surfer/Grimson_new_phylogeny/data/other_data/reformatted_anc_dictU_shuffled'

gene_dict=cPickle.load(open(unshuffled_file, 'rb'))
output_dict=copy.deepcopy(gene_dict)

def shuffle_coords(size):
	split_coords=[[],[]]
	if size%2==0:
		half_point=int(size/2)
	else:
		half_point=int(size/2)+1
#	print size, half_point, len(split_coords[0]), len(split_coords[1])
	for coord in range(0, size, 5):
		index=int(coord/half_point)
		split_coords[index].append([coord, min(coord+5, size)])
	if abs(len(split_coords[0])-len(split_coords[1]))>1:
		print 'error', size, half_point, len(split_coords[0]), len(split_coords[1])
		exit()
	new_coords=[]
	if len(split_coords[1])>len(split_coords[0]):
		print len(split_coords[1]), len(split_coords[0])
	for coord_number, first_half in enumerate(split_coords[0]):
		new_coords.append(first_half)
		if coord_number<len(split_coords[1]):
			new_coords.append(split_coords[1][coord_number])
	return new_coords
	
for gene in gene_dict:
	for ancestor_number, ancestor in enumerate(gene_dict[gene]):
		for descendant_number, descendant in enumerate(gene_dict[gene][ancestor]):
			if ancestor_number==0 and descendant_number==0:
				coord_lists=shuffle_coords(len(gene_dict[gene][ancestor][descendant][0]))
			output_dict[gene][ancestor][descendant][0]=''
			output_dict[gene][ancestor][descendant][1]=''
			for coord in coord_lists:
				start, end=coord
				output_dict[gene][ancestor][descendant][0]+=gene_dict[gene][ancestor][descendant][0][start:end]
				output_dict[gene][ancestor][descendant][1]+=gene_dict[gene][ancestor][descendant][1][start:end]
cPickle.dump(output_dict, open(shuffled_file, 'wb'))
