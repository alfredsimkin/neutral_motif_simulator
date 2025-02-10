'''
a library of modules I use frequently in computational biology. Unlike custom,
this version is modified for use with python3
'''

def count_in_base(string, base, valuestring):
	outlist=[]
	full_length=len(string)
	for count in range(base**len(string)):
		outlist.append(string)
		if string==valuestring[-2]*full_length:
			break
		position=len(string)-1
		value=valuestring.index(string[position])
		string=string[:position]+valuestring[value+1]+string[position+1:]
		while value>=(base-1) and count<(base**len(string)):
			string=string[:position]+valuestring[0]+string[position+1:]
			position-=1
			value=valuestring.index(string[position])
			string=string[:position]+valuestring[value+1]+string[position+1:]
	return outlist

#def revcom(sequence, sequence_type='DNA'):
#	'''
#	the guts of this module are not of my creation
#	I've added some tweaks
#	'''
#	import string
#	if sequence_type=='RNA':
#		complement = string.maketrans('ACGTUN', 'UGCAAN')
#	elif sequence_type=='DNA':
#		complement = string.maketrans('ACGTUN', 'TGCAAN')
#	return sequence.upper().translate(complement)[::-1]

def revcom(seq,nuc='DNA'):
	if nuc=='DNA':
		complement={'A':'T','C':'G','G':'C','T':'A'}
	else:
		complement={'A':'U','C':'G','G':'C','U':'A'}
	return ''.join(reversed([complement[base] for base in seq]))

def meanvar(value_list, vartype='sample'):
	import math
	n=len(value_list)
	total=math.fsum(value_list)
	square_sum=math.fsum([value**2 for value in value_list])
	mean=total/n
	second=2*mean*total
	third=mean**2*n
	if n>1 and vartype=='sample':
		variance=(square_sum-(2*mean*total)+(mean**2)*n)/(n-1)
	else:
		variance=(square_sum/n)-(mean**2)

	return [mean, variance]

def grabcolumn(path, column):
	outlist=[]
	for line in open(path):
		split_line=line.strip().split()
		outlist.append(split_line[column])
	return outlist

def gather_columns(path, columns=None, split_thing='\t'):
	newlist=[]
	for line in open(path):
		line_split=line.strip().split(split_thing)
		if columns==None:
			newline=line_split
		else:
			newline=[line_split[column] for column in columns]
		newlist.append(newline)
	return newlist

def read_fasta(fasta_file):
	seq, name_list, seq_list, seq_dict='', [], [], {}
	for line in open(fasta_file):
		line=line.strip()
		if '>' in line:
			name_list.append(line[1:])
			if len(seq)>0:
				seq_list.append(seq)
				seq=''
		else:
			seq=seq+line
	seq_list.append(seq)
#	for seq_number, name in enumerate(name_list):
#		seq_dict[name]=seq_list[seq_number]
	return [[name, seq_list[name_number]] for name_number, name in enumerate(name_list)]
def print_fasta(fasta_list, outfile, mode='w', line_chars=60):
	"""
	this program prints fasta paired lists to fasta format
	"""
	output=open(outfile, mode)
	for sequence in fasta_list:
		output.write(">"+sequence[0]+"\n")
		for char_number, char in enumerate(sequence[1]):
			output.write(char)
			if char_number%line_chars==line_chars-1 or char_number==(len(sequence[1]))-1:
				output.write("\n")
def print_phylip(species_names, sequence_list, phylipfile, length=500):
	#prints paired name/sequence pairs in sequence_list to an interleaved output file
	#(phylipfile) and names each with names from species_names. Unlike the paired
	#fasta_list, the paired sequence_list has len(species_names) sequences per
	#pairing instead of one
	output=open(phylipfile, 'w')
	total_length=len(sequence_list[1][0])
	output.write(str(len(species_names))+' '+str(total_length)+'\n')
	for species_number, species in enumerate(species_names):
		output.write(species_names[species_number]+' '+sequence_list[1][species_number][:length]+'\n')
	output.write('\n')
	for letter_number, letter in enumerate(sequence_list[1][0][length:]):
		if letter_number%length==0:
			for sequence in sequence_list[1]:
				sequence=sequence[length:]
				output.write(sequence[letter_number:letter_number+length]+'\n')
			output.write('\n')
def print_alignment(fasta_thing, block_size=60):
	#reads an aligned fasta file and prints it 
	#interleaved in blocks of size block_size
	if type(fasta_thing) is str:
		fasta_thing=read_fasta(fasta_thing)
	else:
		fasta_list=fasta_thing
	for sequence in fasta_list:
		print(sequence[0])
	for letter_number, letter in enumerate(fasta_list[0][1]):
		if letter_number%block_size==0:
			for sequence in fasta_list:
				print(sequence[1][letter_number:letter_number+block_size])
			print("\n",)

def split_to_nodes (nodes, unsplit_path, output_folder, master_path):
	import subprocess
	import math
	output_number=1
	new_master_command=open(master_path, 'w')
	unsplit_file=open(unsplit_path)
	job_list=unsplit_file.readlines()
	total_jobs=float(len(job_list))
	divisor=math.ceil(total_jobs/nodes)
	for number, line in enumerate(job_list):
		if number%divisor==0:
			outfile=open("{0}{1}.sh".format(output_folder, output_number), 'w')
			new_master_command.write("{0}{1}.sh\n".format(output_folder, output_number))
			output_number+=1
		if number%divisor==divisor-1 or number==total_jobs-1:
			subprocess.call(["chmod", "+x", "{0}{1}.sh".format(output_folder, output_number-1)])
		outfile.write(line)
	return output_number-1
	
def theta_w(individuals, segsites, sites_w_info):
	#this program takes the number of individuals in the sample,
	#the number of segregating sites in the sample,
	#and the number of total sites considered to calculate theta_w
	a_1=0
	S=float(segsites)/sites_w_info #bases theta on total sites considered rather than total sites in locus
	for value in range(1,individuals):
		a_1=a_1+1.0/value
	theta=S/a_1
	return theta

def replace_all(input_list, conversion_list):
	#this program replaces all first elements from
	#conversion_list found in input_list with all
	#second elements from conversion_list
	for converting_line in conversion_list:
		for in_number, in_line in enumerate(input_list):
			if converting_line[0] in in_line:
				input_list[in_number]=in_line.replace(converting_line[0], converting_line[1])
	return input_list

def overlapping_search(search, string):
	#this definition searches 'string' using search term 'search'
	#and returns the coordinates of all overlapping matches, allowing for gaps
	import re
	expression, starts, ends="", [], []
	gap="-*"
	for letter in search[:-1]:
		expression+=letter+gap
	expression+=search[-1]
	for m in re.finditer("(?="+expression+")", string):
		start=m.start()
		starts.append(start)
		ends.append((re.match(expression, string[start:])).end()+start)
	return [str(start)+':'+str(ends[start_number]) for start_number, start in enumerate(starts)]

def old_choose(n, r):
	#counts ways of choosing r things from a set of n
	import math
	from operator import mul
	if n<r or n<0 or r<0:
		return 0
	if n-r>r:
		factorial_divisor=r
	else:
		factorial_divisor=n-r
	num_list=range(1, n+1)[:-(factorial_divisor+1):-1]
	denominator=math.factorial(factorial_divisor)
	if r>0 and r<n:
		numerator=reduce(mul, num_list)
	else:
		numerator=denominator
	return numerator/denominator

def choose(n, k):
    """
    A fast way to calculate binomial coefficients by Andrew Dalke, from StackOverflow.
    """
    if 0 <= k <= n:
        ntok = 1
        ktok = 1
        for t in xrange(1, min(k, n - k) + 1):
            ntok *= n
            ktok *= t
            n -= 1
        return ntok // ktok
    else:
        return 0
def binomial_prob(prob, number_observed, number_draws):
	'''
	calculates the binomial probability of seeing number_observed events in 
	number_draws total events when the probability of a single event is prob
	'''
	from decimal import Decimal
	ways=Decimal(choose(number_draws, number_observed))
	return Decimal(str(prob**number_observed*(1-prob)**(number_draws-number_observed)))*Decimal(ways)

def gcf(smaller, bigger):
	while smaller>0:
		smaller, bigger=bigger%smaller, smaller
	return bigger

def binomial_prob2(prob_num, prob_den, number_observed, number_draws, fraction=False):
	'''
	a modified (hopefully more precise) means of calculating binomial probabilities
	'''
	from decimal import Decimal
	
	ways=choose(number_draws, number_observed)
	not_prob_num=prob_den-prob_num
	not_observed=number_draws-number_observed
	prob_num=prob_num**number_observed*not_prob_num**not_observed*ways
	prob_den=prob_den**number_observed*prob_den**not_observed
	if fraction:
		return [prob_num, prob_den]
	else:
		return Decimal(prob_num)/Decimal(prob_den)

def cum_binom_prob(prob_num, prob_den, number_observed, number_draws):
	'''
	returns the probability of observing number_observed or more successes in
	number_draws trials given success probability prob_num/prob_den
	'''
	from decimal import Decimal
	if number_observed==0:
		return Decimal(1)
	else:
		for number in range(number_observed):
			if number==0:
				num, den=binomial_prob2(prob_num, prob_den, number, number_draws, 1)
			else:
				new_num, new_den=binomial_prob2(prob_num, prob_den, number, number_draws, 1)
				num, den=num*new_den+new_num*den, den*new_den
				factor=gcf(num, den)
				num, den=num/factor, den/factor
		return Decimal(den-num)/Decimal(den)
	
def check_overlap(start1, end1, start2, end2):
	#checks whether the start and end from the first two entries
	#overlaps with the start and end from the second two entries
	#returns "status" of overlap and start and end coordinates of overlapping bases
	if start1<=start2:
		bigger_small=start2
		status='smaller'
	else:
		bigger_small=start1
		status='bigger'
	if end1<=end2:
		smaller_big=end1
	else:
		smaller_big=end2
	if smaller_big>bigger_small:
		status='overlapping'
	return bigger_small, smaller_big, status

def serpentine(ranked_list, divisor):
	#takes a ranked list and distributes it to divisor
	#sub lists of approximately equal average ranks
	out_list, number=[], 0
	for element in range(divisor):
		out_list.append([])
	for rank_number, value in enumerate(ranked_list):
		if rank_number%divisor==0:
			step=0
		elif (rank_number/divisor)%2==0:
			step=1
		elif (rank_number/divisor)%2==1:
			step=-1
		number=number+step
		out_list[number].append(value)
	return out_list
def fasta_to_fastq(fasta_file_name):
	#converts fasta files into fastq, only works if 
	#each fasta entry is on a single line. All fastq 
	#entries are given a quality score of 'z'
	fastq_file=open(fasta_file_name+'.fq', 'w')
	for line_number, line in enumerate(open(fasta_file_name)):
		if line_number%2==0 and '>' in line:
			fastq_file.write('@'+line[1:])
		elif line_number%2==0 and '>' not in line:
			print('bad fasta file for fastq conversion')
			exit()
		elif line_number%2==1:
			fastq_file.write(line)
			fastq_file.write('+\n')
			fastq_file.write('z'*(len(line)-1)+'\n')
	fastq_file.close()
def fastq_to_fasta(fastq_file):
	#converts fastq files to fasta
	fasta_file=open(fastq_file+'.fa', 'w')
	for line_number, line in enumerate(open(fastq_file)):
		if line_number%4==0 and '@' in line:
			fasta_file.write('>'+line[1:])
		elif line_number%4==0 and '@' not in line:
			print('error')
			exit()
		elif line_number%4==1:
			fasta_file.write(line)
def locate_differences(file1, file2):
	'''
	finds first line where two files differ
	'''
	file1=open(file1)
	file2=open(file2)
	file1_line=file1.readline()
	file2_line=file2.readline()
	line_number=0
	while file1_line==file2_line:
		file1_line=file1.readline()
		file2_line=file2.readline()
		line_number+=1
	print(file1_line)
	print(file2_line)
	print(line_number)
def compute_pi(sequences):
	#computes the popgen statistic pi for a list of DNA sequences
	good_nucs='ACGTacgt'
	diff_count, bad_count=0, 0
	for sequence_number, sequence_one in enumerate(sequences):
		for sequence_two in sequences[sequence_number+1:]:
			for bp_number, bp_one in enumerate(sequence_one):
				bp_two=sequence_two[bp_number]
				if bp_two!=bp_one and bp_one in good_nucs and bp_two in good_nucs:
					diff_count+=1
				if bp_one not in good_nucs or bp_two not in good_nucs:
					bad_count+=1
	total_comparisons=float(choose(len(sequences), 2))*len(sequences[0])
	if bad_count<total_comparisons:
		return diff_count/(total_comparisons-bad_count), total_comparisons-bad_count
	else:
		return -0.02, total_comparisons-bad_count
def compute_pi2(sequences):
	#computes the popgen statistic pi for a list of DNA sequences
	good_nucs='ACGTacgt'
	diff_count, bad_count=0, 0
	for sequence_number, sequence_one in enumerate(sequences):
		for sequence_two in sequences[sequence_number+1:]:
			for bp_number, bp_one in enumerate(sequence_one):
				bp_two=sequence_two[bp_number]
				if bp_two!=bp_one and bp_one in good_nucs and bp_two in good_nucs:
					diff_count+=1
				if bp_one not in good_nucs or bp_two not in good_nucs:
					bad_count+=1
	total_comparisons=float(choose(len(sequences), 2))*len(sequences[0])
	return diff_count, total_comparisons, bad_count
def translate(genetic_sequence):
	#translates some DNA or RNA sequence into peptides
	translated=''
	translations={'UUU':'F', 'CUU':'L', 'AUU':'I', 'GUU':'V', 'UUC':'F', 'CUC':'L',
	'AUC':'I', 'GUC':'V', 'UUA':'L', 'CUA':'L', 'AUA':'I', 'GUA':'V', 'UUG':'L',
	'CUG':'L', 'AUG':'M', 'GUG':'V', 'UCU':'S', 'CCU':'P', 'ACU':'T', 'GCU':'A',
	'UCC':'S', 'CCC':'P', 'ACC':'T', 'GCC':'A', 'UCA':'S', 'CCA':'P', 'ACA':'T',
	'UCG':'S', 'CCG':'P', 'ACG':'T', 'GCG':'A', 'UAU':'Y', 'CAU':'H', 'AAU':'N',
	'GAU':'D', 'UAC':'Y', 'CAC':'H', 'AAC':'N', 'GAC':'D', 'UAA':'*', 'CAA':'Q',
	'AAA':'K', 'GAA':'E', 'UAG':'*', 'CAG':'Q', 'AAG':'K', 'GAG':'E', 'UGU':'C',
	'CGU':'R', 'AGU':'S', 'GGU':'G', 'UGC':'C', 'CGC':'R', 'AGC':'S', 'GGC':'G',
	'UGA':'*', 'CGA':'R', 'AGA':'R', 'GGA':'G', 'UGG':'W', 'CGG':'R', 'AGG':'R',
	'GGG':'G', 'GCA':'A'}
	if 'T' in genetic_sequence:
		genetic_sequence=genetic_sequence.replace('T', 'U')
	for start in range(len(genetic_sequence)//3):
		codon=genetic_sequence[start*3:start*3+3]
		translated+=translations[codon]
	return translated
def make_ranks(x, y):
	'''
	takes paired x and y values and returns ranks, with correct tie handling
	'''
	import scipy.stats
	x=scipy.stats.stats.rankdata(x)
	y=scipy.stats.stats.rankdata(y)
	return x, y

def round_scientific(input_number, precision=3):
	'''
	rounds extremely small or large numbers in scientific notation intelligently
	'''
	number_str=str(input_number)
	number_list=number_str.split('e')
	number=number_list[0]
	if len(number_list)>1:
		exp=number_list[1]
	number=str(round(float(number), precision))
	if len(number_list)>1:
		return float('e'.join([number, exp]))
	else:
		return float(number)

def scatter_it(first, second, first_label, second_label, lims=False, faded=True):
	'''
	takes a list of x values ('first') and y values ('second'), along with x
	labels and y labels, and plots them as a scatter plot. Unless specified,
	limits are auto generated, and data points are faded. Figure will be saved
	as x_label_vs_y_label.png. A linear regression line is also plotted and
	listed as the title of the graph.
	'''
	import scipy.stats
	import matplotlib as mpl
	mpl.use('Agg') #comment this in if using a nongraphical interface, comment out to enable graphical interface
	import matplotlib.pyplot as plt
	import numpy
	m, b, r, linear_p_value, std_err=scipy.stats.linregress(first, second)
#	r, rank_p_value=scipy.stats.spearmanr(first, second)
	y=[m*x+b for x in first]
	if lims!=False:
		plt.xlim(lims[0])
		plt.ylim(lims[1])
	if faded:
		plt.scatter(first, second, color='red', linewidth=0, alpha=0.08)
	else:
		plt.scatter(first, second, color='red', linewidth=0)		
	plt.plot(first, y, color="black")
#	plt.tick_params(labelsize=17)
	plt.xlabel(first_label)
	plt.ylabel(second_label)
#	ax = plt.gca()
#	ax.tick_params(direction='out')
#	plt.xlabel(first_label, size=15)
#	plt.ylabel(second_label, size=15)
#	linear_p_value=numpy.round(linear_p_value, )
	plt.title('m='+str(round_scientific(m))+' b='+str(round_scientific(b))+' r2='+str(round_scientific(r**2))+' p='+str(round_scientific(linear_p_value)), size=15)
#	plt.title('m='+str(round_scientific(m))+' b='+str(round_scientific(b))+' r2='+str(round_scientific(r**2))+' p='+str(round_scientific(rank_p_value)), size=15)
#	plt.grid(True, linestyle='-')
	plt.savefig(first_label+'_vs_'+second_label+'.png', bbox_inches='tight')
	plt.clf()

def MannU(x_values, y_values):
	'''
	takes some x_values and y_values, ranks them together, and outputs which
	has a lower sum of ranks, and a P-value
	'''
	import scipy.stats
	merged_values=x_values+y_values	
	ranks, junk=make_ranks(merged_values, y_values)
	x_ranks=ranks[:len(x_values)]
	rank_sum=sum(x_ranks)
	rank_count=len(x_ranks)
	x_U=rank_sum-((rank_count*(rank_count+1))/2)
	y_U=len(x_values)*len(y_values)-x_U
	smaller_U, p=scipy.stats.mannwhitneyu(x_values, y_values)
	if x_U<y_U:
		statement='x_smaller'
	if y_U<x_U:
		statement='y_smaller'
	if y_U==x_U:
		statement='equal'
	return statement, smaller_U, p

def hypergeometric_prob(pop_success, pop_size, number_observed, number_draws, fraction=False):
	'''
	returns the probability of observing exactly number_observed successes in
	number_draws trials given a total of pop_success successes and a total
	population of pop_size
	'''
	from decimal import Decimal
	success_ways=choose(pop_success, number_observed)
	fail_ways=choose(pop_size-pop_success, number_draws-number_observed)
	numerator=success_ways*fail_ways
	total_ways=choose(pop_size, number_draws)
	if fraction==False:
		return Decimal(numerator)/Decimal(total_ways)
	else:
		return numerator, total_ways

def cum_hypergeometric_prob(pop_success, pop_size, number_observed, number_draws, extra=False):
	from decimal import Decimal
	'''
	returns the probability of observing number_observed successes or more in
	number_draws. Results should be	identical to a one-sided Fisher's exact
	test. (or a two sided if extra=True)
	'''
	num, den=0,1
	for number in range(number_observed):
		new_num, new_den=hypergeometric_prob(pop_success, pop_size, number, number_draws, 1)
		num, den=num*new_den+new_num*den, den*new_den
		factor=gcf(num, den)
		num, den=num/factor, den/factor
	regular=Decimal(den-num)/Decimal(den)
	if extra:
		enrich_value=Decimal(den-num)/Decimal(den)
		exact_num, exact_den=hypergeometric_prob(pop_success, pop_size, number_observed, number_draws, 1)
		deplete_value=Decimal(num)/Decimal(den)+Decimal(exact_num)/Decimal(exact_den)
		if enrich_value<deplete_value:
			return enrich_value, 'enriched'
		else:
			return deplete_value, 'depleted'
	else:
		return regular

def make_float(input_list):
	'''
	converts all numbers in a list of strings (like those read in from input
	files) into floats. Like map(float, input_list) but no error if input isn't
	float
	'''
	for column_number, column in enumerate(input_list):
		try:
			column=float(column)
			input_list[column_number]=column
		except ValueError:
			pass
	return input_list

def get_median(input_list):
	middle_value=len(input_list)/2
	if len(input_list)%2==0:
		median=(input_list[middle_value]+input_list[middle_value-1])/2.0
	else:
		median=input_list[middle_value]
	return median

def parse_csv_line(input_line):
	parsed_list, protected=[],False
	input_line=input_line.strip()
	while len(input_line)>0:
		for character_number, character in enumerate(input_line):
			if character=='"':
				if protected==False:
					protected=True
				else:
					protected=False
					parsed_list.append(input_line[1:character_number])
					input_line=input_line[character_number+2:]
					break
			if character==',' and protected==False:
				parsed_list.append(input_line[:character_number])
				input_line=input_line[character_number+1:]
				break
		if ',' not in input_line and '"' not in input_line:
			parsed_list.append(input_line)
			input_line=''
	return parsed_list

def product(input_list):
	product=1
	for item in input_list:
		product*=item
	return product

def iterate(nested_input_list):
	'''
	takes an arbitrarily deeply nested input list and iterates through all of
	the items in it
	'''
	sizes, final_list=[],[[],[]]
	for small_list in nested_input_list:
		sizes.append(len(small_list))
	total_size=product(sizes)
	for number in xrange(total_size):
		index_list,value_list=[],[]
		for size_number, current_size in enumerate(sizes):
			correction=product(sizes[size_number+1:])
			index=number/correction%current_size
			index_list.append(index)
			value_list.append(nested_input_list[size_number][index])
		final_list[0].append(index_list)
		final_list[1].append(value_list)
	return final_list

