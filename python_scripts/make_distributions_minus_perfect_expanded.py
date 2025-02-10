'''
compares how often every 8mer changes in reality to a distribution of change
rates in simulations. Outputs the number of standard deviations above or below
(+, -) the mean expected simulated rate that observed events occur.
'''
print 'consolidating replicates, this could take some time...'
import custom
import cPickle
import subprocess
import sys
full_dict={}
data_folder, replicate_start, replicate_end, mask_file=sys.argv[1:]
#data_folder=outfolder.split('/')[:-2]
#data_folder='/'.join(data_folder)+'/'
infolder=data_folder+'../simulations/'
mask_dict=cPickle.load(open(mask_file, 'rb'))
real_counts=cPickle.load(open(data_folder+'real_counts'))

def modify_counts(eightmer_counts, anc_seq, desc_seq, gene):
	anc_seq=anc_seq.upper()
	desc_seq=desc_seq.upper()
	for bp_number, bp in enumerate(anc_seq):
		anc_8mer=anc_seq[bp_number:bp_number+8]
		eightmer_set=set(range(bp_number, bp_number+8))
		if len(eightmer_set&mask_dict[gene])==0:
#			print 'mask without overlap'
			continue
#		print 'good'
		if len(anc_8mer)==8 and anc_8mer in good_set:
			desc_8mer=desc_seq[bp_number:bp_number+8]
			if desc_8mer in good_set:
				if anc_8mer not in eightmer_counts:
					eightmer_counts[anc_8mer]=[0,0,0]
				if desc_8mer not in eightmer_counts:
					eightmer_counts[desc_8mer]=[0,0,0]
				eightmer_counts[desc_8mer][0]+=1
				if desc_8mer!=anc_8mer:
					eightmer_counts[anc_8mer][2]+=1
					eightmer_counts[desc_8mer][1]+=1
	return eightmer_counts
eightmer_list=custom.count_in_base('AAAAAAAA', 4, 'ACGTz')
good_set=set(eightmer_list)

for eightmer in real_counts:
	real_counts[eightmer].append(real_counts[eightmer][1]+real_counts[eightmer][2])
for replicate in range(int(replicate_start), int(replicate_end)):
	print 'rep is', replicate
	gene_dict=cPickle.load(open(infolder+'formatted_sim_wrong_anc_dict'+str(replicate), 'rb'))
	sorted_genes=sorted(gene_dict.keys())
	eightmer_counts={}
	for gene_number, gene in enumerate(sorted_genes):
		if gene_number%100==0:
			print gene_number/float(len(sorted_genes))
		for ancestor in gene_dict[gene]:
			for descendant in gene_dict[gene][ancestor]:
				anc, desc=gene_dict[gene][ancestor][descendant]
				eightmer_counts=modify_counts(eightmer_counts, anc, desc, gene)
	for eightmer in eightmer_counts:
		if 'N' in eightmer:
			print eightmer, 'is weird!!! stop!'
		if eightmer not in full_dict:
			full_dict[eightmer]=[[],[],[],[]]
		for column_number, column in enumerate(eightmer_counts[eightmer]):
			full_dict[eightmer][column_number].append(eightmer_counts[eightmer][column_number])
		full_dict[eightmer][3].append(eightmer_counts[eightmer][1]+eightmer_counts[eightmer][2])
subprocess.call(['mkdir', data_folder])
cPickle.dump(full_dict, open(data_folder+'full_sim_counts', 'wb'))
#exit()
real={}
full_count=0
mean_counts={}
print 'done, starting next'
for eightmer in full_dict:
	for column_number, column_values in enumerate(full_dict[eightmer]):
		if column_number not in real:
			real[column_number]={}
			mean_counts[column_number]={}
		mean, variance=custom.meanvar(column_values)
		mean_counts[column_number][eightmer]=mean
		std_dev=variance**0.5
		if eightmer in real_counts:
			if std_dev>0:
				real[column_number][eightmer]=(real_counts[eightmer][column_number]-mean)/std_dev
			elif real_counts[eightmer][column_number]==mean:
				real[column_number][eightmer]=0.0
			elif real_counts[eightmer][column_number]>mean:
				real[column_number][eightmer]='inf'
			else:
				real[column_number][eightmer]='-inf'
real_total=open(data_folder+'real_total', 'w')
real_gains=open(data_folder+'real_gain', 'w')
real_losses=open(data_folder+'real_loss', 'w')
real_sum=open(data_folder+'real_sum', 'w')
for eightmer in real[0]:
	real_total.write(eightmer+'\t'+str(real[0][eightmer])+'\t'+str(real_counts[eightmer][0])+'\t'+str(mean_counts[0][eightmer])+'\n')
for eightmer in real[1]:
	real_gains.write(eightmer+'\t'+str(real[1][eightmer])+'\t'+str(real_counts[eightmer][1])+'\t'+str(mean_counts[1][eightmer])+'\n')
for eightmer in real[2]:
	real_losses.write(eightmer+'\t'+str(real[2][eightmer])+'\t'+str(real_counts[eightmer][2])+'\t'+str(mean_counts[2][eightmer])+'\n')
for eightmer in real[3]:
	real_sum.write(eightmer+'\t'+str(real[3][eightmer])+'\t'+str(real_counts[eightmer][3])+'\t'+str(mean_counts[3][eightmer])+'\n')
