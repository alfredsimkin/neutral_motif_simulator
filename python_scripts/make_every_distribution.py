'''
compares how often every 8mer changes in reality to a distribution of change
rates in simulations. Outputs the number of standard deviations above or below
(+, -) the mean expected simulated rate that observed events occur. Retrofitted
to do this for all seed comparisons.
'''
print 'consolidating replicates, this could take some time...'
import custom
import cPickle
import subprocess
import sys
full_dict={}
data_folder, replicate_start, replicate_end, seed_comparison=sys.argv[1:]
infolder=data_folder+'../simulations/'
outfolder=data_folder
real_counts=cPickle.load(open(data_folder+'everything_real_counts'))
for eightmer in real_counts[seed_comparison]:
	real_counts[seed_comparison][eightmer].append(real_counts[seed_comparison][eightmer][1]+real_counts[seed_comparison][eightmer][2])
for replicate in range(int(replicate_start), int(replicate_end)):
	print 'rep is', replicate
	eightmer_dict=cPickle.load(open(infolder+'everything_simulated_counts'+str(replicate)))
	for eightmer in eightmer_dict[seed_comparison]:
		if 'N' in eightmer:
			print eightmer, 'is weird!!! stop!'
		if eightmer not in full_dict:
			full_dict[eightmer]=[[],[],[],[]]
		for column_number, column in enumerate(eightmer_dict[seed_comparison][eightmer]):
			full_dict[eightmer][column_number].append(eightmer_dict[seed_comparison][eightmer][column_number])
		full_dict[eightmer][3].append(eightmer_dict[seed_comparison][eightmer][1]+eightmer_dict[seed_comparison][eightmer][2])
cPickle.dump(full_dict, open(outfolder+'full_sim_counts_'+seed_comparison, 'wb'))
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
		if eightmer in real_counts[seed_comparison]:
			if std_dev>0:
				real[column_number][eightmer]=(real_counts[seed_comparison][eightmer][column_number]-mean)/std_dev
			elif real_counts[seed_comparison][eightmer][column_number]==mean:
				real[column_number][eightmer]=0.0
			elif real_counts[seed_comparison][eightmer][column_number]>mean:
				real[column_number][eightmer]='inf'
			else:
				real[column_number][eightmer]='-inf'
real_total=open(outfolder+'real_total'+'_'+seed_comparison, 'w')
real_gains=open(outfolder+'real_gain'+'_'+seed_comparison, 'w')
real_losses=open(outfolder+'real_loss'+'_'+seed_comparison, 'w')
real_sum=open(outfolder+'real_sum'+'_'+seed_comparison, 'w')
for eightmer in real[0]:
	real_total.write(eightmer+'\t'+str(real[0][eightmer])+'\t'+str(real_counts[seed_comparison][eightmer][0])+'\t'+str(mean_counts[0][eightmer])+'\n')
for eightmer in real[1]:
	real_gains.write(eightmer+'\t'+str(real[1][eightmer])+'\t'+str(real_counts[seed_comparison][eightmer][1])+'\t'+str(mean_counts[1][eightmer])+'\n')
for eightmer in real[2]:
	real_losses.write(eightmer+'\t'+str(real[2][eightmer])+'\t'+str(real_counts[seed_comparison][eightmer][2])+'\t'+str(mean_counts[2][eightmer])+'\n')
for eightmer in real[3]:
	real_sum.write(eightmer+'\t'+str(real[3][eightmer])+'\t'+str(real_counts[seed_comparison][eightmer][3])+'\t'+str(mean_counts[3][eightmer])+'\n')
