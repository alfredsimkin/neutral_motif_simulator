'''
Attempts to convert unix line endings to DOS line endings, for windows users.
Also adds file extensions, for easier opening on windows.
'''
import sys
import os
import subprocess
folder=sys.argv[1]
subprocess.call(['mkdir', folder+'dos_version/'])

def unix_to_dos(input_file_name, output_file_name):
	'''
	converts unix line endings of an input plain text file to dos line endings
	of an output plain text file
	'''
	output_file=open(output_file_name, 'w')
	for line in open(input_file_name):
		line=line.replace('\n', '\r\n')
		output_file.write(line)
	output_file.close()

for filename in os.listdir(folder):
	if '.' not in filename:
		if 'distribution_' in filename or 'mir_ranks' in filename or 'real' in filename or 'gains_vs_losses' in filename:
			unix_to_dos(folder+filename, folder+'dos_version/'+filename+'.tsv')
		if 'p_values' in filename:
			unix_to_dos(folder+filename, folder+'dos_version/'+filename+'.txt')

