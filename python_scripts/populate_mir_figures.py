'''
First, counts turnover events for all miRNAs (in all species and all
replicates) and sends to a dictionary. Next, averages across replicates and
species to get an 'expected' turnover value for each gene, and across species
in the real dataset to get an 'observed' value. Finally, subtracts expected
from observed, and populates the coordinates from make_empty_mir_figure.py 
'''
import cPickle
import copy
import custom
import sys
import subprocess
data_folder, simulation_folder, real_alignment=sys.argv[1:]
figure_folder=data_folder+'figures/'
subprocess.call(['mkdir', figure_folder])
empty_figures=cPickle.load(open(data_folder+'empty_figures', 'rb'))
real_8mers=set(custom.count_in_base('AAAAAAAA', 4, 'ACGTz'))
real_mirs=set([mir[0] for mir in empty_figures.keys()])

def modify_counts(eightmer_counts, anc_seq, desc_seq):
	anc_seq=anc_seq.upper()
	desc_seq=desc_seq.upper()
	for bp_number, bp in enumerate(anc_seq):
		anc_8mer=anc_seq[bp_number:bp_number+8]
		if len(anc_8mer)==8 and anc_8mer in real_8mers:
			desc_8mer=desc_seq[bp_number:bp_number+8]
			if desc_8mer in real_8mers and (anc_8mer in real_mirs or desc_8mer in real_mirs):
				if anc_8mer not in eightmer_counts:
					eightmer_counts[anc_8mer]=[0,0,0]
				if desc_8mer not in eightmer_counts:
					eightmer_counts[desc_8mer]=[0,0,0]
				eightmer_counts[desc_8mer][0]+=1
				if desc_8mer!=anc_8mer:
					eightmer_counts[anc_8mer][2]+=1
					eightmer_counts[desc_8mer][1]+=1
	return eightmer_counts

def make_full_dict(alignment_file, name, full_dict):
	gene_dict=cPickle.load(open(alignment_file, 'rb'))
	for gene_number, gene in enumerate(gene_dict):
		if gene_number%1000==0:
			print gene_number/float(len(gene_dict))
		for ancestor in gene_dict[gene]:
			for descendant in gene_dict[gene][ancestor]:
				anc, desc=gene_dict[gene][ancestor][descendant]
				real_counts=modify_counts({}, anc, desc)
				for eightmer in real_counts:
					if gene not in full_dict:
						full_dict[gene]={}
					if name not in full_dict[gene]:
						full_dict[gene][name]={}
					if ancestor not in full_dict[gene][name]:
						full_dict[gene][name][ancestor]={}
					if descendant not in full_dict[gene][name][ancestor]:
						full_dict[gene][name][ancestor][descendant]={}
					full_dict[gene][name][ancestor][descendant][eightmer]=copy.deepcopy(real_counts[eightmer])
	return full_dict

def make_all_full_dicts():
	full_dict={}
	print 'getting real'
	full_dict=make_full_dict(real_alignment, 'real', full_dict)
	for sim_number in xrange(100):
		print 'getting sim', sim_number
		full_dict=make_full_dict(simulation_folder+'formatted_sim_wrong_anc_dict'+str(sim_number), 'sim'+str(sim_number), full_dict)
	cPickle.dump(full_dict, open(data_folder+'full_counts', 'wb'), -1)
	return full_dict

def make_species_sum(input_dict):
	summed_dict={}
	species_count=0
	for ancestor in input_dict:
		for descendant in input_dict[ancestor]:
			species_count+=1
			for eightmer in input_dict[ancestor][descendant]:
				if eightmer not in summed_dict:
					summed_dict[eightmer]=[0,0,0]
				for value_number, value in enumerate(input_dict[ancestor][descendant][eightmer]):
					summed_dict[eightmer][value_number]+=value
	return summed_dict, species_count

def make_rep_sum(current_sum, new_sums):
	for eightmer in new_sums:
		if eightmer not in current_sum:
			current_sum[eightmer]=[0,0,0]
		for value_number, value in enumerate(new_sums[eightmer]):
			current_sum[eightmer][value_number]+=value
	return current_sum

def make_average(summed_dict, count):
	average_dict={}
	for eightmer in summed_dict:
		avg=[0,0,0]
		for value_number, value in enumerate(summed_dict[eightmer]):
			avg[value_number]=value/float(count)
		average_dict[eightmer]=avg
	return average_dict

def average_values(full_dict):
	real_avg, sim_avg={},{}
	for gene in full_dict:
#		print gene
		if 'real' in full_dict[gene]:
			real_summed_dict, real_species_count=make_species_sum(full_dict[gene]['real'])
			real_avg[gene]=make_average(real_summed_dict, real_species_count)
		rep_sum={}
		for rep in xrange(100):
			if 'sim'+str(rep) in full_dict[gene]:
				sim_species_sum, junk=make_species_sum(full_dict[gene]['sim'+str(rep)])
				rep_sum=make_rep_sum(rep_sum, sim_species_sum)
		sim_avg[gene]=make_average(rep_sum, real_species_count*100)		
	cPickle.dump(real_avg, open(data_folder+'real_avg_turnover', 'wb'), -1)
	cPickle.dump(sim_avg, open(data_folder+'sim_avg_turnover', 'wb'), -1)	
	return real_avg, sim_avg

def test_presence(input_dict, gene, eightmer):
	result=[0,0,0]
	if gene in input_dict and eightmer[0] in input_dict[gene]:
		result=input_dict[gene][eightmer[0]]
	return result

def populate_figures(real_avg, sim_avg):
#	gene_dict=cPickle.load(open(real_alignment, 'rb'))
	full_figures={}
	for eightmer in empty_figures:
		full_figures[eightmer]=copy.deepcopy(empty_figures[eightmer])
		for gene in empty_figures[eightmer][1]:
			coords=empty_figures[eightmer][1][gene]
#			print eightmer, gene, coords
			real=test_presence(real_avg, gene, eightmer)
			sim=test_presence(sim_avg, gene, eightmer)
			full_figures[eightmer][1][gene]={}
			full_figures[eightmer][1][gene]['gain']=coords+[real[1]-sim[1]]
			full_figures[eightmer][1][gene]['loss']=coords+[real[2]-sim[2]]
#			if gene=='NM_005578' and eightmer[0]=='GTGCCATA':
#				print eightmer[0], gene, 'gain', real[1]-sim[1], real[1], sim[1]
#				exit()
#			if real[1]-sim[1]!=0 and eightmer[0]=='ACCAAAGA':
#				print eightmer[0], gene, 'gain', real[1]-sim[1], real[1], sim[1]
#			if real[2]-sim[2]!=0 and eightmer[0]=='ACCAAAGA':
#				print eightmer[0], gene, 'loss', real[2]-sim[2], real[2], sim[2]
	return full_figures

make_all_full_dicts()
full_dict=cPickle.load(open(data_folder+'full_counts', 'rb'))
real_avg, sim_avg=average_values(full_dict)
#real_avg=cPickle.load(open(data_folder+'real_avg_turnover', 'rb'))
#sim_avg=cPickle.load(open(data_folder+'sim_avg_turnover', 'rb'))
full_figures=populate_figures(real_avg, sim_avg)
cPickle.dump(full_figures, open(figure_folder+'populated_figures', 'wb'), -1)
