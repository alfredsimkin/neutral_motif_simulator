'''
takes an input file of introns, keeps introns belonging to genes analyzed in 3'
UTRs, and, because introns are so much larger than 3' UTRs, prints coords of
some fraction of the center of each intron using a scaling factor (100=1/100th)
Output is sent to bed_files for find_maf_overlaps to work on.

sep. 2019 update: to make this fully compatible with all_together.py,
attempting to create a bed file of introns that can then be run through the
full all_together.py script
'''

import cPickle
import subprocess
data_folder='/mnt/hulk_tank/alfred_stuff/silver_surfer/Grimson_new_phylogeny/data/other_data/9-25-19_introns/'
intron_file='/home/alfred/refseq_hg38_introns.bed'
good_dict=cPickle.load(open('/mnt/hulk_tank/alfred_stuff/silver_surfer/Grimson_new_phylogeny/data/other_data/reformatted_anc_dictU', 'rb'))
strand_bed='/mnt/hulk_tank/alfred_stuff/silver_surfer/Grimson_new_phylogeny/data/other_data/refseq_hg38.bed'
long_beds=[]
scale_factor=10
good_genes=set(good_dict.keys())
found_genes=set([])
found_coords=set([])
outdict={}

def make_strand_dict(strand_bed):
	'''
	records the strand that every gene falls on and outputs as a dictionary
	'''
	strand_dict={}
	for line in open(strand_bed):
		line=line.strip().split()
		gene='_'.join(line[3].split('_')[:2])
		if gene in good_genes:
			strand=line[5]
			if gene not in strand_dict:
				strand_dict[gene]=strand
			elif gene in strand_dict and strand!=strand_dict[gene]:
				print line, 'is weird!!!'
				exit()
	return strand_dict

def get_hundredth(line, gene):
	chrom, start, end=line[0], int(line[1]), int(line[2])
	length=(end-start)/scale_factor
	middle=(end+start)/2
	line[1]=middle-length
	line[2]=middle+length
	line[3]=gene
	if line[2]>line[1]:
		return '\t'.join(map(str, line))
	else:
		return []

def make_long_bed(short_line, strand):
	'''
	modification that creates fake bed entries with parameters that should
	allow the intron entries to be counted as 3' UTR material having a single
	exon (for use in all_together.py).
	'''
	line=short_line.split('\t')
	chrom, start, end, gene=line[0], int(line[1]), int(line[2]), line[3]
	rgb, score, exons, exon_starts=0, 0, 1, '0,'
	exon_lengths=str(end-start)+','
	if strand=='+':
		new_start=start-1
		new_end=end
		coding_start=new_start
		coding_end=start
	elif strand=='-':
		new_start=start
		new_end=end+1
		coding_start=end
		coding_end=new_end
	else:
		print strand
	long_bed='\t'.join(map(str, [chrom, new_start, new_end, gene, score, strand, coding_start, coding_end, rgb, exons, exon_lengths, exon_starts]))
	return long_bed

strand_dict=make_strand_dict(strand_bed)

for org_line in open(intron_file):
	split_line=org_line.strip().split()
	gene='_'.join(split_line[3].split('_')[:2])
	chrom, start, end=split_line[:3]
	if chrom not in outdict and 'hap' not in chrom:
		outdict[chrom]=[]
	coords=(chrom, start, end)
	if gene in good_genes and coords not in found_coords:
		line=get_hundredth(split_line, split_line[3])
		if line and 'hap' not in chrom:
			outdict[chrom].append(line)
			found_coords.add(coords)
			found_genes.add(gene)
			strand=strand_dict[gene]
			long_beds.append(make_long_bed(line, strand))
#commented out lines below can all be commented back in as a block to restore original functionality of this program
#subprocess.call(['mkdir', '-p', data_folder+'bed_files/'])
#output2=open(data_folder+'bed_files/chrom_list', 'w')
#output3=open(data_folder+'revised_introns.bed', 'w')
output4=open(data_folder+'revised_introns_as_3UTRs.bed', 'w')
output4.write('\n'.join(long_beds))
#for chrom in outdict:
#	if 'hap' not in chrom:
#		output_file=open(data_folder+'bed_files/'+chrom+'.bed', 'w')
#		output2.write(chrom+'\n')
#		for gene in outdict[chrom]:
#			output_file.write(gene+'\n')
#			output3.write(gene+'\n')
#		output_file.close()
