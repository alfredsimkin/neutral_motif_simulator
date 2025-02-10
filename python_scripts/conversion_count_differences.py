'''
calculates the difference between the number of observed sites and the number of expected sites for every miRNA seed
'''

def print_table(key_list, conversion_dict):
	print '\t'.join(key_list)
	for number in xrange(75):
		line_list=[]
		for key in key_list:
			line_list+=[str(sorted(conversion_dict[key])[number])]
		print '\t'.join(line_list)

data_folder='/home/alfred/Dropbox/Simkin_Grimson_final_ms_folder/figures_and_data/mammal_data/deeply_conserved_mirs(human_mouse_fish)/conversions_0_100_none_6mer_7mer1A_7merm8_8mer_human_mouse_fish_best_testing_pipeline/'
site_types=['none', '6mer', '7mer1A', '7merm8', '8mer']
conversion_dict={}
gain_keys, loss_keys=[],[]
for number_one, first_site in enumerate(site_types):
	for second_site in site_types[number_one+1:]:
		gain_key=first_site+'->'+second_site
		loss_key=second_site+'->'+first_site
		gain_keys+=[gain_key]
		loss_keys+=[loss_key]
		for line in open(data_folder+'mir_ranks_gain_'+first_site+'_'+second_site):
			line=line.strip().split('\t')
			expected_gains=float(line[-2])
			observed_gains=float(line[-3])
			gain_diff=[observed_gains-expected_gains]
			if gain_key not in conversion_dict:
				conversion_dict[gain_key]=[]
			conversion_dict[gain_key]+=gain_diff
		for line in open(data_folder+'mir_ranks_loss_'+first_site+'_'+second_site):
			line=line.strip().split('\t')
			expected_losses=float(line[-2])
			observed_losses=float(line[-3])
			loss_diff=[observed_losses-expected_losses]
			if loss_key not in conversion_dict:
				conversion_dict[loss_key]=[]
			conversion_dict[loss_key]+=loss_diff
print 'gains'
print_table(gain_keys, conversion_dict)
print 'losses'
print_table(loss_keys, conversion_dict)
