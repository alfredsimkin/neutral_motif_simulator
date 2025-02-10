'''
really a special use program designed to remove redundant miRNAs from multiple
runs of 'populate_real_mers.py.'
'''

included_files=['human_mirs', 'mouse_mirs']
excluded_files=['chicken_mirs', 'fish_mirs']
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

for mir in mir_dict['human_mirs']:
	if mir in intersected_mirs:
		output_file.write(mir_dict['human_mirs'][mir])
