'''
compares how often every 8mer changes in reality to a distribution of change
rates in simulations. Outputs the number of standard deviations above or below
(+, -) the mean expected simulated rate that observed events occur.

this version (expanded) adds counts to the real_total, real_gain, and real_loss output files
'''
print 'consolidating replicates, this could take some time...'
import custom
import cPickle
import subprocess
import sys
full_dict={}
data_folder, replicate_start, replicate_end=sys.argv[1:]
#data_folder=outfolder.split('/')[:-2]
#data_folder='/'.join(data_folder)+'/'
infolder=data_folder
outfolder=data_folder
real_counts=cPickle.load(open(data_folder+'real_branch_counts'))

for ancestor in real_counts:
	for descendant in real_counts[ancestor]:
		for eightmer in real_counts[ancestor][descendant]:
			real_counts[ancestor][descendant][eightmer].append(real_counts[ancestor][descendant][eightmer][1]+real_counts[ancestor][descendant][eightmer][2]) #append a 'sum' count to the existing gain and loss counts

for replicate in range(int(replicate_start), int(replicate_end)):
	print 'rep is', replicate
	eightmer_dict=cPickle.load(open(infolder+'sim_branch_counts'+str(replicate)))
	for ancestor in eightmer_dict:
		if ancestor not in full_dict:
			full_dict[ancestor]={}
		for descendant in eightmer_dict[ancestor]:
			if descendant not in full_dict[ancestor]:
				full_dict[ancestor][descendant]={}
			print ancestor, descendant
			for eightmer in eightmer_dict[ancestor][descendant]:
				if 'N' in eightmer:
					print eightmer, 'is weird!!! stop!'
				if eightmer not in full_dict[ancestor][descendant]:
					full_dict[ancestor][descendant][eightmer]=[[],[],[],[]]
				for column_number, column in enumerate(eightmer_dict[ancestor][descendant][eightmer]):
					full_dict[ancestor][descendant][eightmer][column_number].append(eightmer_dict[ancestor][descendant][eightmer][column_number])
				full_dict[ancestor][descendant][eightmer][3].append(eightmer_dict[ancestor][descendant][eightmer][1]+eightmer_dict[ancestor][descendant][eightmer][2])
cPickle.dump(full_dict, open(outfolder+'full_sim_counts', 'wb'))

#exit()
#full_dict=cPickle.load(open(outfolder+'full_sim_counts', 'rb'))

real={} #the badly named 'real' dictionary stores the number of standard deviations above or below the mean of simulated values that the real value falls at
full_count=0
mean_counts={}
print 'done, starting next'
for ancestor in full_dict:
	print ancestor
	if ancestor not in real:
		real[ancestor]={}
		mean_counts[ancestor]={}
	for descendant in full_dict[ancestor]:
		if descendant not in real[ancestor]:
			real[ancestor][descendant]={}
			mean_counts[ancestor][descendant]={}
		for eightmer in full_dict[ancestor][descendant]:
			for column_number, column_values in enumerate(full_dict[ancestor][descendant][eightmer]):
				if column_number not in real[ancestor][descendant]:
					real[ancestor][descendant][column_number]={}
					mean_counts[ancestor][descendant][column_number]={}
				mean, variance=custom.meanvar(column_values)
				mean_counts[ancestor][descendant][column_number][eightmer]=mean
				std_dev=variance**0.5
				if eightmer in real_counts[ancestor][descendant]:
					if std_dev>0:
						real[ancestor][descendant][column_number][eightmer]=(real_counts[ancestor][descendant][eightmer][column_number]-mean)/std_dev
					elif real_counts[ancestor][descendant][eightmer][column_number]==mean:
						real[ancestor][descendant][column_number][eightmer]=0.0
					elif real_counts[ancestor][descendant][eightmer][column_number]>mean:
						real[ancestor][descendant][column_number][eightmer]='inf'
					else:
						real[ancestor][descendant][column_number][eightmer]='-inf'
print 'outputting to output files'
for ancestor in full_dict:
	for descendant in full_dict[ancestor]:
		real_total=open(outfolder+ancestor+'-'+descendant+'real_total', 'w')
		real_gains=open(outfolder+ancestor+'-'+descendant+'real_gain', 'w')
		real_losses=open(outfolder+ancestor+'-'+descendant+'real_loss', 'w')
		real_sum=open(outfolder+ancestor+'-'+descendant+'real_sum', 'w')
		for eightmer in real[ancestor][descendant][0]:
			real_total.write(eightmer+'\t'+str(real[ancestor][descendant][0][eightmer])+'\t'+str(real_counts[ancestor][descendant][eightmer][0])+'\t'+str(mean_counts[ancestor][descendant][0][eightmer])+'\n')
		for eightmer in real[ancestor][descendant][1]:
			real_gains.write(eightmer+'\t'+str(real[ancestor][descendant][1][eightmer])+'\t'+str(real_counts[ancestor][descendant][eightmer][1])+'\t'+str(mean_counts[ancestor][descendant][1][eightmer])+'\n')
		for eightmer in real[ancestor][descendant][2]:
			real_losses.write(eightmer+'\t'+str(real[ancestor][descendant][2][eightmer])+'\t'+str(real_counts[ancestor][descendant][eightmer][2])+'\t'+str(mean_counts[ancestor][descendant][2][eightmer])+'\n')
		for eightmer in real[ancestor][descendant][3]:
			real_sum.write(eightmer+'\t'+str(real[ancestor][descendant][3][eightmer])+'\t'+str(real_counts[ancestor][descendant][eightmer][3])+'\t'+str(mean_counts[ancestor][descendant][3][eightmer])+'\n')
