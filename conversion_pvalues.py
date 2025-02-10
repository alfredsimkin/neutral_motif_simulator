'''
Calculates two types of p-values for site conversions: overall probability of
strengthening compared against overall probability of weakening, and overall
probability of these for each group.
'''

import custom
folders=[['/home/alfred/Dropbox/Simkin_Grimson_final_ms_folder/figures_and_data/mammal_data/deeply_conserved_mirs(human_mouse_fish)/conversions_0_100_none_6mer_7mer1A_7merm8_8mer_human_mouse_fish_best_testing_pipeline/', 1, 'together'], ['/home/alfred/big_data/turnover_rates/4-28-15_deeper_phylogeny/conversions_0_100_none_6mer_7mer1A_7merm8_8mer_human_mouse_fish_best_group1/', 3, 'group1'], ['/home/alfred/big_data/turnover_rates/4-28-15_deeper_phylogeny/conversions_0_100_none_6mer_7mer1A_7merm8_8mer_human_mouse_fish_best_group2/', 3, 'group2'],['/home/alfred/big_data/turnover_rates/4-28-15_deeper_phylogeny/conversions_0_100_none_6mer_7mer1A_7merm8_8mer_human_mouse_fish_best_group3/', 3, 'group3']]
site_types=['none', '6mer', '7mer1A', '7merm8', '8mer']
result_dict={}
for folder in folders:
	gain_list, loss_list=[],[]
	folder_path, column_number, folder_name=folder
	for event in ['mir_ranks_gain_', 'mir_ranks_loss_']:
		for site_number, first_site in enumerate(site_types):
			for second_site in site_types[site_number+1:]:
				file_name=event+first_site+'_'+second_site
				stats_list=[int(line.strip().split('\t')[column_number]) for line in open(folder_path+file_name)]
				if 'gain' in event:
					gain_list.extend(stats_list)
				else:
					loss_list.extend(stats_list)
	gains_vs_losses=custom.MannU(gain_list, loss_list)
	result_dict[folder_name]=[gain_list, loss_list]
