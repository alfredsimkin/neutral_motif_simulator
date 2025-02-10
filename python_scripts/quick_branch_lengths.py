'''
a quick script to retrieve the branch lengths associated with any given list of 8mers
'''
eightmer_set=set(['GCAAAAAA', 'CGAACAAA', 'GCACTTTA', 'GTGCCTTA', 'CTACCTCA', 'TTGCACTA', 'GTGCCAAA', 'TGTTTACA', 'ACCAAAGA', 'GTGCAATA', 'AGCAATAA', 'TGGTGCTA', 'AAGCACAA', 'CTATGCAA', 'AAGCCATA', 'TACTTGAA', 'GACAATCA', 'TTGCCAAA', 'GTGCCATA', 'TTTGCACA', 'CTCAGGGA', 'GGACCAAA', 'GCACCTTA', 'GTACTGTA', 'ACATTCCA', 'TGCTGCTA', 'TTCCGTTA', 'TGAATGTA', 'AATGTGAA', 'ACTTTATA', 'TACGGGTA', 'AGCATTAA', 'GACTGTTA', 'ACTACTGA', 'AGCACTTA', 'CAGTATTA', 'TGCACTGA', 'CGTCTTAA', 'ACATATCA', 'ATGTAGCA', 'ATGCTGCA', 'ATAAGCTA', 'ACTACCTA', 'CACTGTGA', 'ACAGGGTA', 'AAAGGGAA', 'GTCTTCCA', 'AACTGGAA', 'ATGAAGGA', 'CACTGCCA', 'CACCAGCA', 'CAGTGTTA', 'CTGAGCCA', 'CTGTTACA', 'TCATCTCA', 'AGTCTTAA', 'TAGGTCAA', 'AACGGTTA', 'CGCTGCTA', 'ACGCACAA', 'CATTTCAA', 'ACACTCCA', 'AGAGATTA', 'ATGCTGGA', 'GTGTCATA', 'TCCGTCCA', 'AACTGACA', 'TGAGATTA', 'ATGCAGTA', 'CGGTACGA', 'ACTGTAGA', 'AGACACGA', 'CTGTGGTA', 'CCTGCTGA', 'AGTTCTCA'])
#branch_lengths={line.strip().split('\t')[1]:float(line.strip().split('\t')[2]) for line in open('mature_mir_branch_lengths')}

branch_lengths={}
for line in open('mature_mir_branch_lengths'):
	line=line.strip().split('\t')
	if line[1] in eightmer_set:
		if line[1] not in branch_lengths:
			branch_lengths[line[1]]=[]
		branch_lengths[line[1]].append(float(line[2]))

for eightmer in eightmer_set:
	print eightmer, max(branch_lengths[eightmer])
