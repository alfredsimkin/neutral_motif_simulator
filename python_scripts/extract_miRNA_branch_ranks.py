'''
a script to get miRNA ranks from a table containing number of standard
deviations every 8mer falls relative to mean of simulated values
'''
import sys
import custom
data_folder=sys.argv[1]
mirfile=sys.argv[2]
pairing=sys.argv[3]
real_mirs={line.strip().split('\t')[0]:line.strip().split('\t')[1] for line in open(mirfile)}
stat_list=['gain', 'loss', 'sum']

for stat in stat_list:
	rank_list=[]
	output_file=open(data_folder+pairing+'_mir_ranks_'+stat, 'w')
	p_value_file=open(data_folder+pairing+'_p_values_'+stat, 'w')
	all_8mers=[line.strip().split('\t') for line in open(data_folder+pairing+'real_'+stat)]
	all_8mers=sorted([[float(line[1]), line[0], line[2], line[3]] for line in all_8mers])
	for line_number, line in enumerate(all_8mers):
#		print line[1]
		if line[1] in real_mirs:
			output_file.write(line[1]+'\t'+str(line_number)+'\t'+str(line[0])+'\t'+line[2]+'\t'+line[3]+'\t'+real_mirs[line[1]]+'\n')
			rank_list.append(line_number)
	smaller, U, p_value=custom.MannU(rank_list, range(line_number+1))
	if smaller=='x_smaller':
		p_value_file.write(' '.join(['mirs', stat, 'ranks lower than all 8mers,', 'p-value:', str(p_value)])+'\n')
	if smaller=='y_smaller':
		p_value_file.write(' '.join(['mirs', stat, 'ranks higher than all 8mers,', 'p-value:', str(p_value)])+'\n')
	rank_list=map(str, rank_list)
	p_value_file.write('ranks listed below:\n')
	p_value_file.write('\n'.join(rank_list))
