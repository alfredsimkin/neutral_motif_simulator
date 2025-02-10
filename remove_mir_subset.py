'''
really a special use program designed to combine miRNAs from multiple
runs of 'populate_real_mers2.py' and to remove redundant miRNAs from other
runs of this program.

This version (2) uses the power of sets to specifically subtract particular
miRNA groups

Dependencies:
1. input files generated by populate_real_mers.py

Input parameters:
1. names of input files from populate_real_mers.py whose miRNA seeds should be
included in the final list if present in all of these input files.
2. names of input files from populate_real_mers.py whose miRNA seeds should be
excluded from the final list.
'''

#input parameters here:
included_files=['hsamdo_minus_dre_mirs']
excluded_files=['mmu_mirs']

included_name='_'.join([name.split('_')[0] for name in included_files])
excluded_name='_'.join([name.split('_')[0] for name in excluded_files])
if len(excluded_files)>0:
	output_file=open(included_name+'_minus_'+excluded_name+'_mirs', 'w')
else:
	output_file=open(included_name+'_mirs', 'w')

mir_dict={}
for file_number, included_file in enumerate(included_files):
	mir_dict[included_file]={line.strip().split('\t')[0]:line for line in open(included_file)}
	if file_number==0:
		intersected_mirs=set(mir_dict[included_file].keys())
	else:
		intersected_mirs=intersected_mirs&set(mir_dict[included_file].keys())
for excluded_file in excluded_files:
	mir_dict[excluded_file]={line.strip().split('\t')[0]:line for line in open(excluded_file)}
	intersected_mirs=intersected_mirs-set(mir_dict[excluded_file].keys())

for mir in mir_dict[included_files[0]]:
	if mir in intersected_mirs:
		output_file.write(mir_dict[included_files[0]][mir])
