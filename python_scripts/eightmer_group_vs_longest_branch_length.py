'''
gathers the conservation levels of every mature miRNA from
mature_mir_branch_lengths, extracts 8mer seeds, and for each 8mer, reports the
branch length of the mature miRNA that is most conserved. Ranks these 8mers
from least conserved to most conserved, and reports whether the ranks of 8mers
corresponding to miRNAs in group 1 are significantly lower than the ranks of
8mers corresponding to miRNAs in group 2. Also compares group 1 to group 3 and
group 2 to group 3
'''
import custom
root_folder='/home/alfred/Dropbox/Simkin_Grimson_final_ms_folder/figures_and_data/mammal_data/deeply_conserved_mirs(human_mouse_fish)/'
loss_file=root_folder+'main_results(vanilla)/mir_ranks_loss'
gain_file=root_folder+'main_results(vanilla)/mir_ranks_gain'

def get_group(gain_rank, loss_rank):
	if gain_rank<3277 and loss_rank<3277:
		group=1
	elif gain_rank>=3277 and loss_rank<3277:
		group=2
	elif gain_rank>=3277 and loss_rank>=3277:
		group=3
	elif gain_rank<3277 and loss_rank>=3277:
		group=4
	return group

def get_ranks(loss_file, gain_file):
	'''
	extracts ranks from the gain and loss files and returns groups for 8mer
	seeds by calling get_group
	'''
	group_dict={}
	gain_list=[line.strip().split('\t') for line in open(gain_file)]
	loss_list=[line.strip().split('\t') for line in open(loss_file)]
	mirseed, rank, name=0,1,5
	gain_dict={line[mirseed]:[line[rank], line[name]] for line in gain_list}
	loss_dict={line[mirseed]:[line[rank], line[name]] for line in loss_list}
	for mirseed in gain_dict:
		mir_name=gain_dict[mirseed][1]
		gain_rank=int(gain_dict[mirseed][0])
		loss_rank=int(loss_dict[mirseed][0])
		group=get_group(gain_rank, loss_rank)
		group_dict[mirseed]=group
	return group_dict

group_dict=get_ranks(loss_file, gain_file)

cons_dict={}
for line in open('mature_mir_branch_lengths'):
	line=line.strip().split()
	if line[1] not in cons_dict or float(line[2])>cons_dict[line[1]]:
		cons_dict[line[1]]=float(line[2])

group_lengths_dict={}
for mirseed in group_dict:
	group=group_dict[mirseed]
	if group not in group_lengths_dict:
		group_lengths_dict[group]=[]
	group_lengths_dict[group]+=[cons_dict[mirseed]]

group_list=list(range(1,5))
for group_number, first_group in enumerate(group_list):
	for second_group in group_list[group_number+1:]:
		print first_group, second_group
		print custom.MannU(group_lengths_dict[first_group], group_lengths_dict[second_group])
