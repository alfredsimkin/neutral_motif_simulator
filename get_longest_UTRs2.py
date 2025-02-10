#takes the longest transcript of all possible isoforms for a given 3' UTR.
#An older version of this script checks to make sure no 3' UTRs are contained
#entirely within a 3' UTR intron, as these might derive from entire intronic
#genes, but since none of these were found, it seems safe to assume these
#independent, intronic genes don't exist.

#this version additionally has been tweaked to use a bed file with one 3' UTR
#exon per entry as input. Ultimately, the script needs to be merged with a high
#up part of my turnover rate pipeline (one that knows the genomic start
#and end coordinates of every 3' UTR and its strand)

import cPickle
import sys
data_folder, utr_exons=sys.argv[1:]
#data_folder='/media/alfred/data_drive/big_data/turnover_rates/'
#org_gene_dict=cPickle.load(open(data_folder+'reformatted_wrong_anc_dictU', 'rb'))
output_file=open(data_folder+utr_exons[:-4]+'_best_isoform.bed', 'w')
gene_dict, org_gene_dict={},{}

of_interest=['NM_001319139.1', 'NM_001319141.1', 'NM_013411.4', 'NM_001199199.1', 'NM_001319142.1', 'NM_001319140.1', 'NM_001625.3', 'NM_001319143.1']

def sort_bed(bed_name):
	bed_list=[]
	for line in open(bed_name):
		line=line.strip().split()
		rearranged=[line[0]]+[line[5]]+line[1:5]
		newline=[]
		for column in rearranged:
			try:
				column=int(column)
			except ValueError:
				pass
			newline.append(column)
		bed_list.append(newline)
	new_list=sorted(bed_list)
	final_list=[]
	for line in new_list:
		unarranged=[line[0]]+line[2:]+[line[1]]
		final_list.append(unarranged)
	return final_list

def populate_groups(gene_dict):
	group_counter=0
	group_dict={}
	previous_chrom='yellow'
	bed_list=sort_bed(data_folder+utr_exons)
	for line in bed_list:
		gene=line[3]
		chrom, start, end, strand, length=gene_dict[gene]
		coords=tuple(gene_dict[gene][:-1])
		if group_counter>1253 and group_counter<1257:
			print group_counter
			print line
			print chrom, previous_chrom, start, previous_start, end, previous_end, strand, previous_strand
		if chrom!=previous_chrom or start>group_end or strand!=previous_strand:
			if group_counter>1253 and group_counter<1257:
				print 'outside'
			group_counter+=1
			biggest_length=length
			biggest_gene=gene
			group_start=start
			group_end=end
			genes=set([gene])
			group_dict[group_counter]={coords:[genes, [biggest_length, biggest_gene]]}
		elif chrom==previous_chrom and start<group_end and strand==previous_strand:
			if group_counter>1253 and group_counter<1257:
				print 'inside'
			genes=previous_genes|set([gene])
			group_start=min(previous_start, start)
			group_end=max(previous_end, end)
			coords=(chrom, group_start, group_end, strand)
			if length>biggest_length:
				biggest_length=length
				biggest_gene=gene
			group_dict[group_counter]={coords:[genes, [biggest_length, biggest_gene]]}
		previous_chrom, previous_start, previous_end, previous_genes, previous_length, previous_strand=chrom, start, end, genes, length, strand
	return group_dict
		
#populating the gene_dict
for line in open(data_folder+utr_exons):
	org_line=line.strip()
	line=line.strip().split()
	gene=line[3]
	if gene not in org_gene_dict:
		org_gene_dict[gene]=[]
	org_gene_dict[gene].append(org_line)
	if gene not in gene_dict:
		gene_dict[gene]=[line[0], int(line[1]), int(line[2]), line[5], int(line[2])-int(line[1])]
	else:
		gene_dict[gene][4]+=int(line[2])-int(line[1]) #total size of element increases
	if int(line[1])<gene_dict[gene][1]:
		gene_dict[gene][1]=int(line[1])
	if int(line[2])>gene_dict[gene][2]:
		gene_dict[gene][2]=int(line[2])
group_dict=populate_groups(gene_dict)

'''
for group in group_dict:
	for coords in group_dict[group]:
		for gene in group_dict[group][coords][0]:
			if gene in of_interest:
				print group
				print group_dict[group]
exit()

#putting every gene into its own 'group'

group_counter=0
group_dict={}
for gene_number, gene in enumerate(gene_dict):
	if gene_number%100==0:
		print float(gene_number)/len(gene_dict)
	not_found=True
	for group in group_dict:
		group_coords=group_dict[group].keys()[0]
		group_chrom, group_start, group_end, group_strand=group_coords
		if gene_chrom==group_chrom and gene_strand==group_strand and gene_start<group_end and gene_end>group_start:
			not_found=False
			gene_list=group_dict[group][group_coords]
			group_dict[group][group_coords][0].add(gene)
			if length>group_dict[group][group_coords][1][0]:
				group_dict[group][group_coords][1]=[length, gene]
			group_dict[group].pop(group_coords)
			group_dict[group][(gene_chrom, min(group_start, gene_start), max(group_end, gene_end), group_strand)]=gene_list
	gene_chrom, gene_start, gene_end, gene_strand, length=gene_dict[gene]
	gene_coords=tuple(gene_dict[gene][:-1])
	group_dict[group_counter]={gene_coords:[set([gene]), [length, gene]]}
	group_counter+=1
'''
chrom, start, end, strand=0,1,2,3
new_size=2
old_size=12
groups=group_dict.keys()
group_start=0
while old_size!=new_size:
	groups=groups[group_start:]
	print old_size, new_size
	old_size=len(group_dict)
	new_size=old_size #test line, I think this should fix an infinite run loop in situations that don't need collapsing to the longest isoform without breaking situations that do need collapsing
	broken=False
	for group_number, group1 in enumerate(groups):
		group_start=group_number
		if group_number%100==0:
			print group_number/float(len(groups))
		for group2 in groups[group_number+1:]:
			one_coords=group_dict[group1].keys()[0]
			two_coords=group_dict[group2].keys()[0]
			one_set=group_dict[group1][one_coords][0]
			two_set=group_dict[group2][two_coords][0]
			if len(one_set&two_set)>0:
				print group_dict[group1][one_coords]
				print group_dict[group1][one_coords][1][0]
				if group_dict[group1][one_coords][1][0]>group_dict[group2][two_coords][1][0]:
					new_biggest=group_dict[group1][one_coords][1]
				else:
					new_biggest=group_dict[group2][two_coords][1]
				group_dict.pop(group2)
				groups.remove(group2)
				new_coords=(one_coords[chrom], min(one_coords[start], two_coords[start]), max(one_coords[end], two_coords[end]), one_coords[strand])
				group_dict[group1].pop(one_coords)
				group_dict[group1][new_coords]=[one_set|two_set, new_biggest]
				new_size=len(group_dict)
				broken=True
				break
		if broken:
			break

print 'final size was', new_size
new_spliced={}
for group in group_dict:
	coords=group_dict[group].keys()[0]
	best_transcript=group_dict[group][coords][1][1]
	for line in org_gene_dict[best_transcript]:
		output_file.write(line+'\n')
