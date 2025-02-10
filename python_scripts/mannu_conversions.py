import custom
states=['none','6mer','7mer1A','7merm8','8mer']
data_folder='/mnt/hulk_tank/alfred_stuff/silver_surfer/Grimson_new_phylogeny/data/other_data/conversions_0_100_none_6mer_7mer1A_7merm8_8mer_human_mouse_fish_group1/'
output_file=open('mannu_conversions', 'w')
for first_state_number, first_state in enumerate(states):
	for second_state in states[first_state_number+1:]:
		loss_ranks=[int(line.strip().split()[1]) for line in open(data_folder+'mir_ranks_loss_'+first_state+'_'+second_state)]
		gain_ranks=[int(line.strip().split()[1]) for line in open(data_folder+'mir_ranks_gain_'+first_state+'_'+second_state)]
		state, u, p=custom.MannU(loss_ranks, range(65536))
		stringed=map(str, [state, u, p, len(loss_ranks)*65536, first_state, second_state, 'loss'])
		output_file.write('\t'.join(stringed)+'\n')
		state, u, p=custom.MannU(gain_ranks, range(65536))
		stringed=map(str, [state, u, p, len(gain_ranks)*65536, first_state, second_state, 'gain'])
		output_file.write('\t'.join(stringed)+'\n')
