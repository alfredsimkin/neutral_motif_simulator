#unlike previous versions of this script, this version outputs all 3' UTR
#exons, does not output intronic material, and does not deal with fasta files.
#this version fixes a bug that caused redundantly named transcripts with
#multiple genomic locations (only one of which had a 5' UTR) to be missed by my
#filter that removes redundantly located genes. The bug affects only 4 genes in
#5' UTRs and is irrelevant for 3' UTRs as it doesn't affect the output.

import sys
import custom
data_folder, refseq_input_bed=sys.argv[1:]
BED_gene_list=[line.strip().split('\t') for line in open(data_folder+refseq_input_bed)]
def extract_boundaries(gene):
	gene_length=0
	genomic_start, genomic_end=gene[1:3]
	strand, coding_start, coding_end=gene[5:8]
	has_UTR, coding_coords, UTR_5_coords, UTR_3_coords=0, [], [], []
	gene_name, gene_chrom=gene[3], gene[0]
	if coding_end>coding_start:
		if strand=='-' and genomic_start<coding_start:
			has_UTR='yes'
		if strand=='+' and genomic_end>coding_end:
			has_UTR='yes'
	if has_UTR and 'hap' not in gene_chrom and coding_end>coding_start:
		starts_list=[int(start)+genomic_start for start in gene[11].split(',')[:-1]]
		lengths_list=[int(length) for length in gene[10].split(',')[:-1]]
		for start_number, start in enumerate(starts_list):
			end=start+lengths_list[start_number]
			if end>coding_start and start<coding_start:
				if strand=='-':
					UTR_3_coords.append([start, coding_start])
			if end>coding_end and start<coding_end:
				if strand=='+':
					UTR_3_coords.append([coding_end, end])
			if (end<=coding_start and strand=='-') or (start>=coding_end and strand=='+'):
				UTR_3_coords.append([start, end])
			gene_length+=(end-start)
	return [gene_name, gene_chrom, strand, UTR_3_coords, gene_length]

def convert_to_int(str_list):
	for string_number, string in enumerate(str_list):
		try:
			str_list[string_number]=int(str_list[string_number])
		except ValueError:
			pass
	return str_list

unique_set, non_unique_set, converted_list, good_genes=set([]), set([]), [], set([])

for gene in BED_gene_list:
	gene=convert_to_int(gene)
	gene_attributes=extract_boundaries(gene)
#	print gene_attributes
	gene_name, gene_chrom, UTR_coords=gene_attributes[0], gene_attributes[1], gene_attributes[3]
	if UTR_coords:
		converted_list.append(gene_attributes)
	if 'hap' not in gene_chrom:
		if gene_name not in unique_set:
			unique_set.add(gene_name)
		else:
			non_unique_set.add(gene_name)

UTR_introns=open(data_folder+refseq_input_bed.replace('.bed', '')+'_3UTR_introns.bed', 'w')
UTR_positions=open(data_folder+refseq_input_bed.replace('.bed', '')+'_3UTR_exons.bed', 'w')
for converted_gene in converted_list:
	gene_name, gene_chrom, strand, UTR_3, gene_length=converted_gene
	if gene_name not in non_unique_set and 'hap' not in gene_chrom:
		for UTR_exon in UTR_3:
			UTR_positions.write('\t'.join(map(str, [gene_chrom, UTR_exon[0], UTR_exon[1], gene_name, 0, strand]))+'\n')
		if len(UTR_3)>=1:# and len(fasta_gene_dict[gene_name])==gene_length:
			good_genes.add(gene_name)
			if len(UTR_3)>1:
				entry_number=1
				for entry in UTR_3[1:]:
					UTR_introns.write('\t'.join([gene_chrom, str(UTR_3[entry_number-1][1]), str(entry[0]), gene_name])+'\n')
					entry_number+=1
print len(good_genes)
filtered_bed=open(data_folder+refseq_input_bed.replace('.bed', '')+'_filtered.bed', 'w')
for bed in BED_gene_list:
	if bed[3] in good_genes:
		filtered_bed.write('\t'.join(map(str, bed))+'\n')
