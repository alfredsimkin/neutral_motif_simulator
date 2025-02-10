'''
A script to analyze the effects of single gene removals. Calculates how many
genes have any effect when removed, how many have a positive effect, how many
have a negative effect, and the cumulative effect of gene removal when summed

unlike analyze_single_removals2.py, this version prints the effect of each gene
removal in terms of difference from the mean of simulations rather than change
in standard deviations

This version compares each replicate against all other replicates, as well as
the real against all replicates but the first (to have identical n)

This version (2) plots standard deviations instead of differences
'''
import custom
import sys
import math
import cPickle
import copy
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

input_folder, mirfile, data_folder=sys.argv[1:]
output_folder=input_folder
stat_types=['gains', 'losses']
all_eightmers=set([line.strip().split('\t')[0] for line in open(mirfile)])
all_genes=cPickle.load(open(data_folder+'all_gene_names', 'rb'))

def fill_zeroes(input_dict): #may be redundant
	for gene in all_genes:
		for eightmer in all_eightmers:
			if gene not in input_dict:
				input_dict[gene]={}
			if eightmer not in input_dict[gene]:
				input_dict[gene][eightmer]=[0,0]
	return input_dict

def populate_master_dict():
	full_sim_counts={}
	for rep in range(0, 100): #hard coded, should be passed
		rep_sim_counts=cPickle.load(open(input_folder+'miRNA_sim_counts'+str(rep)))
		rep_sim_counts=fill_zeroes(rep_sim_counts)
		print rep
		for gene in all_genes:
			if gene not in full_sim_counts:
				full_sim_counts[gene]={}
			for eightmer in all_eightmers:
				if eightmer not in full_sim_counts[gene]:
					full_sim_counts[gene][eightmer]=[[],[]]
				if gene in rep_sim_counts and eightmer in rep_sim_counts[gene]:
					full_sim_counts[gene][eightmer][0].append(rep_sim_counts[gene][eightmer][0])
					full_sim_counts[gene][eightmer][1].append(rep_sim_counts[gene][eightmer][1])
				else: #probably redundant
					print 'weird'
					full_sim_counts[gene][eightmer][0]=0
					full_sim_counts[gene][eightmer][1]=0
	cPickle.dump(full_sim_counts, open(output_folder+'full_sim_counts', 'wb'), -1)
	return full_sim_counts

def get_median(input_list): #maybe legacy and therefore removable
	size=len(input_list)
	middle=size/2
	if size%2==0:
		median=(input_list[middle]+input_list[middle-1])/2.0
	else:
		median=input_list[middle]
	return median

def get_stats(input_list, eightmer, stat, rep): #maybe legacy and therefore removable
	input_list.sort()
	lowest=input_list[0]
	highest=input_list[-1]
	lower=[item for item in input_list if item<0]
	higher=[item for item in input_list if item>0]
	total=len(lower)+len(higher)
	one=total/100
	if total>0 and one==0:
		one=1
	highest_one=higher[-one:]
	lowest_one=lower[:one]
	if len(lower)>0:
		low_mean, low_var=custom.meanvar(lower)
		low_median=get_median(lower)
		lowest_mean, lowest_var=custom.meanvar(lowest_one)
	else:
		low_mean, low_var, low_median, lowest_mean=0,0,0,0
	if len(higher)>0:
		high_mean, high_var=custom.meanvar(higher)
		high_median=get_median(higher)
		highest_mean, highest_var=custom.meanvar(highest_one)
	else:
		high_mean, high_var, high_median, highest_mean=0,0,0,0
	outlist=[eightmer, stat, rep, len(lower), len(higher), low_mean, high_mean, low_median, high_median, lowest, highest, lowest_mean, highest_mean]
	print '\t'.join(map(str, outlist))

def make_sd(mean, var, score):
	if score==mean:
		sd=0
	elif var==0:
		sd='invalid'
	else:
		sd=(score-mean)/(var**0.5)
	return sd

def get_real_background(real_counts, full_sim_counts):
	output_file=open(output_folder+'real_distributions', 'wb')
	for stat_number, stat in enumerate(stat_types):
		for eightmer in all_eightmers:
			background_list=[]
			for gene in full_sim_counts:
				background=full_sim_counts[gene][eightmer][stat_number][1:]
				foreground=real_counts[gene][eightmer][stat_number]
				mean, var=custom.meanvar(background)
				excess=make_sd(mean, var, foreground)
				if excess!=0 and excess!='invalid':
					background_list.append([excess, gene])
			#outlist=get_stats(background_list, eightmer, stat, 'real')
			#master_out.append(outlist)
			master_out=[eightmer, stat, 'real', background_list]
			cPickle.dump(master_out, output_file) #why is this sequential dumps rather than one big dump?

def get_fake_background(full_sim_counts):
	output_file=open(output_folder+'fake_distributions', 'wb')
	for stat_number, stat in enumerate(stat_types):
		print stat, stat_number
		for eightmer_number, eightmer in enumerate(all_eightmers):
			print eightmer, eightmer_number
			for rep in range(0,100): #hard coded, should be passed
				background_list=[]
				for gene in full_sim_counts:
					background=copy.deepcopy(full_sim_counts[gene][eightmer][stat_number])
					foreground=background.pop(rep)
					mean, var=custom.meanvar(background)
					excess=make_sd(mean, var, foreground)
					if excess!=0 and excess!='invalid':
						background_list.append([excess, gene])
					#outlist=get_stats(background_list, eightmer, stat, rep)
					#master_out.append(outlist)
				master_out=[eightmer, stat, rep, background_list]
				cPickle.dump(master_out, output_file) #why is this sequential dumps rather than one big dump?

full_sim_counts=populate_master_dict()
#full_sim_counts=cPickle.load(open(output_folder+'full_sim_counts', 'rb'))
get_fake_background(full_sim_counts)
real_counts=cPickle.load(open(input_folder+'miRNA_gene_counts'))
real_counts=fill_zeroes(real_counts)
get_real_background(real_counts, full_sim_counts)
'''
cPickle.dump(sim_means, open(output_folder+'sim_count_means', 'wb'), -1)
diff_dict={}
#print all_eightmers
for stat_number, stat_type in enumerate(stat_types):
	diff_dict[stat_type]={}
#	outfile.write('\t'.join(['eightmer', 'slow_genes', 'fast_genes', 'slow_max', 'fast_max', 'slow_mean', 'fast_mean', 'slow_median', 'fast_median'])+'\n')
	for gene in all_genes:
		for eightmer in all_eightmers:
			if eightmer not in diff_dict[stat_type]:
				diff_dict[stat_type][eightmer]={'all':0}
			real_score=real_counts[gene][eightmer][stat_number]
			sim_mean=sim_means[gene][eightmer][stat_number][0]
			sim_stdev=sim_means[gene][eightmer][stat_number][1]**0.5
			if sim_stdev!=0:
				diff=(real_score-sim_mean)/sim_stdev
			elif real_score==sim_mean:
				diff=0
			elif real_score>sim_mean:
				diff=1.5
			elif real_score<sim_mean:
				diff=-1.5
			else:
				print 'weird', real_score, sim_mean, sim_stdev
				exit()
			if diff!=0:
				diff_dict[stat_type][eightmer][gene]=diff
				diff_dict[stat_type][eightmer]['all']+=diff
for stat in diff_dict:
	output_file=open(output_folder+stat+'_statistics', 'w')
	for eightmer in diff_dict[stat]:
		score_list, fast, slow=[],[],[]
		for gene in diff_dict[stat][eightmer]:
			if gene!='all':
				score=diff_dict[stat][eightmer][gene]
				if score<-2:
					slow.append(score)
				elif score>2:
					fast.append(score)
				score_list.append(score)
		score_list.sort()
		if not len(fast)>0:
			fast=[0,0]
		if not len(slow)>0:
			slow=[0,0]
		slow_mean, slow_var=custom.meanvar(slow)
		fast_mean, fast_var=custom.meanvar(fast)
		stringed=map(str, [eightmer, diff_dict[stat][eightmer]['all'], len(slow), len(fast), min(slow), max(fast), slow_mean, fast_mean, get_median(slow), get_median(fast), len(slow)/float(len(fast)), abs(slow_mean)-fast_mean])
		output_file.write('\t'.join(stringed)+'\n')
		plt.title(eightmer)
		plt.plot(range(len(score_list)), score_list)
		plt.plot(range(len(score_list)), [0 for number in range(len(score_list))])
		plt.margins(0,0)
		plt.savefig(output_folder+eightmer+'_'+stat)#, bbox_inches='tight')
		plt.clf()
'''
