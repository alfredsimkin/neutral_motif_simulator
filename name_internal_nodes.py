import sys
data_folder, input_tree=sys.argv[1:]
output_file=open(data_folder+input_tree.replace('.nh', '')+'mod.nh', 'w')

input_string=''.join([line.strip() for line in open(data_folder+input_tree)])
input_list=input_string.split('):')
new_string=''
for item_number, item in enumerate(input_list[:-1]):
	new_string+=item+')'+str(item_number+1)+':'
new_string+=input_list[-1]
output_file.write(new_string+'\n')
