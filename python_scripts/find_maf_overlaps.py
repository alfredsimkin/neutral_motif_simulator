'''
input files:
bed folder containing:
bed formatted files of coordinate ranges, one file per chrom
chrom_list (list of chromosomes)
maf folder containing:
maf formatted files, one file per chrom
'''
import sys
data_folder, maf_folder, species=sys.argv[1:]
good_species=set(eval(species))
def get_next_two(byte_in):
	two_file.seek(byte_in)
	maf_entry=[]
	maf_line=two_file.readline()
	while maf_line:
		if len(maf_line)>1 and maf_line.startswith('s '):
			split_maf_line=maf_line.split()
			if split_maf_line[1].split('.')[0] in good_species:
				maf_entry.append(split_maf_line[1:])
		elif len(maf_line)==1:
			if len(maf_entry)==0:
				print 'maf file breaks here at byte count',
				print two_file.tell()
				exit()
			break
		maf_line=two_file.readline()
	byte_count=two_file.tell()
	if maf_entry:
		chrom=maf_entry[0][0].split('.')[1]
		name='place_holder'
		maf_start=int(maf_entry[0][1])
		maf_end=maf_start+int(maf_entry[0][2])
		attributes=[chrom, name, maf_start, maf_end]
	else:
		attributes=[]
	return attributes, byte_count, maf_entry

def get_next_one(byte_in):
	one_file.seek(byte_in)
	line=one_file.readline().strip().split()
	if len(line)>1:
		line[1], line[2]=int(line[1]), int(line[2])
		line=[line[0], line[3], line[1], line[2]]
	byte_count=one_file.tell()
	return line, byte_count

def analyze_overlappers(two_values, one_values, chunk):
	human_seq=chunk[0][5]
	desired=(max(one_values[2], two_values[2]), min(one_values[3], two_values[3]))
	current_pos=one_values[2]
	started, ended='no','no'
	for pos, letter in enumerate(human_seq):
		if started=='no' and current_pos>=desired[0]:
			started=pos
		if ended=='no' and current_pos>=desired[1]:
			ended=pos
			break
		if letter!='-':
			current_pos+=1
	if started=='no' and current_pos>=desired[0]:
		started=pos+1
	if ended=='no' and current_pos>=desired[1]:
		ended=pos+1
	overlapper=[]
	for sequence in chunk:
		overlapper.append([sequence[0], sequence[5][started:ended]])
	return desired, overlapper

def create_indices(get_next, filetype):
	start, indices, line=0, [], 'not_empty'
	counter=0
	while line!=[]:
		indices.append(start)
		line, start=get_next(start)[:2]
#		print start
		if counter%100000==0:
			print indices[-1]/1000000, 'megabytes of', filetype, 'file indexed so far'
		counter+=1
	return indices[:-1]

import cPickle
bed_folder=data_folder+'bed_files/'
chroms=open(bed_folder+'chrom_list')
#maf_folder='/media/alfred/data_drive/big_data/joint_grad_work/multiz_100way/'
#modify these two files to match your pre-sorted inputs
utr_dict={}
dict_counter, bed_counter=0,0
for chromosome in chroms:
	chromosome=chromosome.strip()
	print '\n\n\n', chromosome
	reset_counter='None'
	try:
		two_file=open(maf_folder+chromosome+'.maf')
	except IOError:
		print '*********\n\n'+maf_folder+chromosome+'.maf', 'does not seem to exist, skipping\n\n**********'
		continue
	print 'making indices, reading MAF file'
	two_indices=create_indices(get_next_two, 'MAF')
#	cPickle.dump(one_indices, open('one_indices', 'w'))
#	one_indices=cPickle.load(open('one_indices'))
	print '\nmaf file read, reading BED file'
	one_file=open(bed_folder+chromosome+'.bed')
	one_indices=create_indices(get_next_one, 'BED')
	print 'bed file read'
	print 'there are', len(one_indices), 'bed entries and', len(two_indices), 'maf chunks to be analyzed'
	one_count, two_count=0,0
	print 'comparing indexed files to extract coordinates of interest...'
	while 1:
#		if one_count>1640 and one_count<1660 and chromosome=='chr11':
#			interesting=True
#		else:
		interesting=False
		if two_count<len(two_indices):
			two_values, byte_count, chunk=get_next_two(two_indices[two_count])
			two_chrom, two_name, two_start, two_end=two_values
		if one_count<len(one_indices):
			one_values=get_next_one(one_indices[one_count])[0]
			one_chrom, one_name, one_start, one_end=one_values
		if one_chrom<two_chrom or (one_chrom==two_chrom and one_end<two_start) or two_count>len(two_indices)-1: #one is smaller, look at next one, rewind two_count
			if interesting:
				print '\none is smaller'
				print one_name, one_start, one_end
				print two_name, two_start, two_end
			if one_count>=len(one_indices)-1:
				break
			if reset_counter!='None':
				two_count=reset_counter
			bed_counter+=1
			if bed_counter%3000==2999:
				print 'dumping to', dict_counter
				cPickle.dump(utr_dict, open(data_folder+'chunk_dict'+str(dict_counter), 'wb'), -1)
				dict_counter+=1
				utr_dict={}
			if one_count%100==0:
				print 'fraction of', chromosome, 'analyzed is', one_count/float(len(one_indices))
			one_count+=1
		elif two_chrom<one_chrom or (one_chrom==two_chrom and two_end<one_start): #two is smaller, advance two, advance rewind_counter
			if interesting:
				print '\none is bigger'
				print one_name, one_start, one_end
				print two_name, two_start, two_end
			two_count+=1
			if reset_counter!='None':
				reset_counter+=1
		elif one_chrom==two_chrom and one_end>=two_start and one_start<=two_end:
			if interesting:
				print '\noverlapping'
				print one_name, one_start, one_end
				print two_name, two_start, two_end
			if reset_counter=='None': #they overlap, analyze and initialize a rewind counter
				reset_counter=two_count
			overlap_tuple, overlap_chunk=analyze_overlappers(one_values, two_values, chunk)
			if one_name not in utr_dict:
				utr_dict[one_name]={}
			utr_dict[one_name][overlap_tuple]=overlap_chunk
			two_count+=1
		if two_count>=len(two_indices)-1 and one_count>=len(one_indices)-1:
			break
cPickle.dump(utr_dict, open(data_folder+'chunk_dict'+str(dict_counter), 'wb'), -1)
max_file=open(data_folder+'max_rep', 'w')
max_file.write(str(dict_counter)+'\n')
