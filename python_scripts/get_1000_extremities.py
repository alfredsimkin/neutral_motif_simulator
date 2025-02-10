'''
Gets the 1000 most extreme 8mers for producing figure 7

Run 'remove_known_motifs6.py' on the output to produce data for figure 7
'''

gain_file='/mnt/hulk_tank/alfred_stuff/silver_surfer/Grimson_new_phylogeny/data/other_data/paper_stats/real_gain'
loss_file='/mnt/hulk_tank/alfred_stuff/silver_surfer/Grimson_new_phylogeny/data/other_data/paper_stats/real_loss'
top_loss_file=open('top_1000_losses', 'w')
top_gain_file=open('top_1000_gains', 'w')
bottom_loss_file=open('bottom_1000_losses', 'w')
bottom_gain_file=open('bottom_1000_gains', 'w')

def get_fraction(start, end, input_file, output_file):
	new_list=[line.strip().split()[:2] for line in open(input_file)]
	new_list=sorted([[float(line[1]), line[0]] for line in new_list])
	sliced=new_list[start:end]
	for item in sliced:
		item=[item[1], item[0]]
		output_file.write('\t'.join(map(str, item))+'\n')

get_fraction(0, 1000, gain_file, top_gain_file)
get_fraction(0, 1000, loss_file, top_loss_file)
get_fraction(-1000, None, gain_file, bottom_gain_file)
get_fraction(-1000, None, loss_file, bottom_loss_file)
