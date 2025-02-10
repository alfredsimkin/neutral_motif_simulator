'''
sums the number of site conversions (expected and observed) for every site conversion type
'''
data_folder='/home/alfred/Dropbox/Simkin_Grimson_final_ms_folder/figures_and_data/mammal_data/deeply_conserved_mirs(human_mouse_fish)/conversions_0_100_none_6mer_7mer1A_7merm8_8mer_human_mouse_fish_best_testing_pipeline/'
site_types=['none', '6mer', '7mer1A', '7merm8', '8mer']
print '\t'.join(['first_site', 'second_site', 'observed_gains', 'expected_gains', 'observed_losses', 'expected_losses'])
for number_one, first_site in enumerate(site_types):
	for second_site in site_types[number_one+1:]:
		expected_gains, observed_gains, expected_losses, observed_losses=0,0,0,0
		for line in open(data_folder+'mir_ranks_gain_'+first_site+'_'+second_site):
			line=line.strip().split('\t')
			expected_gains+=float(line[-2])
			observed_gains+=float(line[-3])
		for line in open(data_folder+'mir_ranks_loss_'+first_site+'_'+second_site):
			line=line.strip().split('\t')
			expected_losses+=float(line[-2])
			observed_losses+=float(line[-3])
		print '\t'.join(map(str, [first_site, second_site, observed_gains, expected_gains, observed_losses, expected_losses]))
