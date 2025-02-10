'''
Takes a user specified miRNA and set of tissues, and gets miRNA turnover rates for each gene in
that tissue (using output from populate_mir_figures.py). Then, takes ranked
gene expression data, both relative and absolute (from make_x_y.py) and graphs
the output.
'''
import cPickle
import custom
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy
from matplotlib import colors

data_folder='/home/alfred/grimson_postdoc/big_data/turnover_rates/coexpression_stuff/old_school/'
mir_turnovers=cPickle.load(open(data_folder+'populated_figures', 'rb'))
mir_exp=cPickle.load(open(data_folder+'miRNA_expression/mir_gene_exp_dict2', 'rb'))
shared_tissue_dict=dict([line.strip().split() for line in open(data_folder+'tissue_key')][1:])
#for key in shared_tissue_dict.keys():
#	value=shared_tissue_dict[key]
#	shared_tissue_dict[value]=key
#	shared_tissue_dict.pop(key)
xy_values=cPickle.load(open(data_folder+'x_and_y_dict', 'rb'))
rel_tissue_exp=cPickle.load(open(data_folder+'human_gene_atlas/'+'gene_exp_dict', 'rb'))
k=0.02

def calculate_distance(first_point, second_point):
	x_distance=abs(first_point[0]-second_point[0])
	y_distance=abs(first_point[1]-second_point[1])
	return (x_distance**2+y_distance**2)**0.5

def distribute_values(current_x, current_y, intensity, value_list):
	for x in range(84):
		for y in range(84):
			r=calculate_distance([current_x, current_y], [x, y])
			result=(1/(r**2+(k*84)))*intensity
			value_list[x][y]+=result
	return value_list

#miRNA=('GTACTGTA', 'hsa-miR-101-3p')
miRNA=('GTGCAATA', 'hsa-miR-32-5p, hsa-miR-92a-3p, hsa-miR-92b-3p, hsa-miR-363-3p, hsa-miR-367-3p')
#miRNA=('CACTGCCA', 'hsa-miR-34a-5p, hsa-miR-449a, hsa-miR-34c-5p, hsa-miR-449b-5p')
#miRNA=('TGGTGCTA', 'hsa-miR-29a-3p, hsa-miR-29b-3p, hsa-miR-29c-3p')
#miRNA=('GTGCCTTA', 'hsa-miR-124-3p, hsa-miR-506-3p')
exp_file=open(miRNA[1]+'_'+'miRNA_rel_exp', 'w')
shared_tissue_list=[]
#tissue_list=xy_values.keys()
#tissue_list=['Testis', 'TestisLeydigCell', 'TestisIntersitial', 'Prostate', 'TestisSeminiferousTubule', 'TestisGermCell']
#tissue='CD34+'
#shared_tissue_list=['Cerebellum', 'Heart', 'PancreaticIslet', 'Uterus', 'PrefrontalCortex', 'Prostate', 'Thyroid', 'Pituitary', 'Testis', 'Placenta', 'Liver', 'CD19+_BCells(neg._sel.)', 'Ovary', 'CD8+_Tcells', 'CD56+_NKCells', 'MedullaOblongata', 'CD4+_Tcells', 'CD14+_Monocytes', 'CD34+']

for tissue in mir_exp[miRNA]:
	print tissue
	if tissue[1] in shared_tissue_dict:
		shared_tissue_list.append(shared_tissue_dict[tissue[1]])
		exp_file.write(str(tissue)+'\n')
exp_file.close()
exit()
for tissue in shared_tissue_list:
	final_values={}
	turnover_dict=mir_turnovers[miRNA][1]
	all_genes=set(turnover_dict.keys())
	for gene_number, gene in enumerate(all_genes):
		if gene_number%1000==0:
			print gene_number
		for event in turnover_dict[gene]:
			if event not in final_values:
				final_values[event]=[[0]*84 for number in range(84)]
			x,y=xy_values[tissue][gene]
			intensity=turnover_dict[gene][event][2]
			final_values[event]=distribute_values(x, y, intensity, final_values[event])
	for event in final_values:
		max_smoothed, min_smoothed=0,0
		for x_number, x in enumerate(final_values[event]):
			for y_number, y in enumerate(x):
				smoothed=y
				if smoothed>max_smoothed:
					max_smoothed=smoothed
				if smoothed<min_smoothed:
					min_smoothed=smoothed
		cmap = colors.ListedColormap([(0.7, 0.1, 0.1, 0.8), (0.7, 0.1, 0.1, 0.6), (0.7, 0.1, 0.1, 0.4), (0.7, 0.1, 0.1, 0.2), (0.7, 0.1, 0.1, 0.1), (0.1, 0.1, 0.7, 0.1), (0.1, 0.1, 0.7, 0.2), (0.1, 0.1, 0.7, 0.4), (0.1, 0.1, 0.7, 0.6), (0.1, 0.1, 0.7, 0.8)])
		smallest=min(abs(max_smoothed), abs(min_smoothed))
		increment=(smallest)/5.0
		bounds=[(-1*smallest)+(value_number*increment) for value_number in range(11)]
		norm=colors.BoundaryNorm(bounds, cmap.N)
		heatmap = plt.pcolor(numpy.array(final_values[event]), cmap=cmap, norm=norm)
		plt.colorbar(heatmap, ticks=bounds)
		plt.savefig(miRNA[1]+'_'+tissue+'_'+event+'.png', bbox_inches='tight')
		plt.clf()
