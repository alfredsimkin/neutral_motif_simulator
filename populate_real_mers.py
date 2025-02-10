'''
Uses a table of reads downloaded from mirbase 21, that has been modified to
label the arms of a miRNA and, where reads are available for both arms,
designate each arm as 'major' or 'minor'. Uses this table and a list of mature
miRNAs to extract seeds that are present (and not the minor isoform) in a fixed
list of species, and prints these seeds along with their names.

Dependencies:
1. custom.py script (provided in this folder)
2. file with read counts for miRNAs (miRNA_read_table.tsv, in this folder)
3. file with mature miRNAs (mirbase_21_mature, in this folder)

Input parameters:
1. names of species to examine
2. name of output file
'''
import custom

#input parameters here:
species_list=['mmu']
output_name='mmu_mirs'

def assign_maj():
	maj_dict={}
	for line in open('miRNA_read_table.tsv'):
		line=line.strip().split()
		if line[4] not in maj_dict:
			maj_dict[line[4]]=[]
		maj_dict[line[4]].append(line[7])
		if line[9] not in maj_dict:
			maj_dict[line[9]]=[]
		maj_dict[line[9]].append(line[12])
	return maj_dict

maj_dict=assign_maj()
mirlist=[line.strip().split() for line in open('mirbase_21_mature')]
output_file=open(output_name, 'w')
mer_dict={}
for line in mirlist:
	mer=custom.revcom(line[-1][1:8])+'A'
	if mer not in mer_dict:
		mer_dict[mer]=[]
	mer_dict[mer].append(line[0])

for mer in mer_dict:
	new_list=[]
	test_dict={species:0 for species in species_list}	
	for mir in mer_dict[mer]:
		if mir[:3] in test_dict:
			if mir[:3]==species_list[0]:
				new_list.append(mir)
			if mir in maj_dict:
				status=maj_dict[mir]
				if len(set(status))>1 or list(set(status))[0]!='minor':
					test_dict[mir[:3]]+=1
	failed=False
	for species in test_dict:
		if test_dict[species]==0:
			failed=True
	if not failed:
		ref_mirs=', '.join(new_list)
		all_mirs=', '.join(mer_dict[mer])
		output_file.write('\t'.join(map(str, [mer, ref_mirs, all_mirs]))+'\n')
