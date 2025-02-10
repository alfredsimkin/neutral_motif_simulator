import subprocess
import sys
data_folder, input_bed=sys.argv[1:]
subprocess.call(['mkdir', data_folder+'bed_files'])

found_coords=set([])
outdict={}

def start_pos(input_line):
	return int(input_line.split('\t')[1])

weird_file=open('weird_file', 'w')
for line in open(data_folder+input_bed):
	split_line=line.strip().split()
	gene=split_line[3]
	chrom, start, end=split_line[:3]
	if chrom not in outdict:
		outdict[chrom]=[]
	coords=(chrom, start, end)
	if coords not in found_coords:
		outdict[chrom].append(line)
		found_coords.add(coords)
	else:
		print 'weird!'
		weird_file.write(line)
#		exit()
for chrom in outdict:
	outdict[chrom]=sorted(outdict[chrom], key=start_pos)
output2=open(data_folder+'bed_files/chrom_list', 'w')
for chrom in outdict:
	output_file=open(data_folder+'bed_files/'+chrom+'.bed', 'w')
	output2.write(chrom+'\n')
	for gene in outdict[chrom]:
		output_file.write(gene)
	output_file.close()
