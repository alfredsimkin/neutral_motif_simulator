'''
creates x and y coordinates for each gene, where x is 'absolute' expression and
y is 'relative' expression
'''
import cPickle
import numpy
import sys
data_folder=sys.argv[1]

probe_gene_dict=cPickle.load(open(data_folder+'gene_exp_dict', 'rb'))
probe_tissue_dict=cPickle.load(open(data_folder+'tissue_exp_dict', 'rb'))

#mir_gene_dict=cPickle.load(open(data_folder+'miRNA_expression/'+'mir_gene_exp_dict', 'rb'))
#mir_tissue_dict=cPickle.load(open(data_folder+'miRNA_expression/'+'mir_tissue_exp_dict', 'rb'))


def get_x(tissue_dict):
	x_coordinates={}
	for tissue in tissue_dict:
		if tissue not in x_coordinates:
			x_coordinates[tissue]={}
		tissue_dict[tissue]
#		print len(tissue_dict)
		divisor=len(tissue_dict[tissue])/float(len(tissue_dict))
		for gene_number, gene in enumerate(tissue_dict[tissue]):
			exp, name=gene
			bin_number=int(gene_number/divisor)
#			print gene, tissue, bin_number
			x_coordinates[tissue][name]=bin_number
	return x_coordinates

def get_y(gene_dict):
	y_coordinates={}
	for gene in gene_dict:
		if gene not in y_coordinates:
			y_coordinates[gene]={}
		for tissue_number, tissue in enumerate(gene_dict[gene]):
			exp, name=tissue
#			print gene, tissue, tissue_number
			y_coordinates[gene][name]=tissue_number
	return y_coordinates

def combine(x_coordinates, y_coordinates):
	final_dict={}
	for tissue in x_coordinates:
		if tissue not in final_dict:
			final_dict[tissue]={}
		for gene in x_coordinates[tissue]:
			x=x_coordinates[tissue][gene]
			y=y_coordinates[gene][tissue]
			final_dict[tissue][gene]=[x,y]
	return final_dict

def analyze_final(final_dict):
	for tissue in final_dict:
		the_dict={}
		for value in range(len(final_dict)):
			counter=0
			for gene in final_dict[tissue]:
				if final_dict[tissue][gene][1]==value:
					counter+=1
			if counter not in the_dict:
				the_dict[counter]=0
			the_dict[counter]+=1
		print the_dict

x_coordinates=get_x(probe_tissue_dict)
y_coordinates=get_y(probe_gene_dict)
final_dict=combine(x_coordinates, y_coordinates)
cPickle.dump(final_dict, open(data_folder+'x_and_y_dict', 'wb'), -1)
print 'done with make_x_y.py'
#analyze_final(final_dict)
#print final_dict.keys()
