import sys
data_folder, good_species, prefix=sys.argv[1:]
good_species=eval(good_species)
species_set=set(good_species)
intron_list=[line.strip().split() for line in open(data_folder+prefix+'_3UTR_introns.bed')]
intron_dict={}
for intron in intron_list:
	if intron[3] not in intron_dict:
		intron_dict[intron[3]]=[]
	intron_dict[intron[3]].append(map(int, intron[1:3]))

import cPickle
gene_tracker={}
special_genes={}
max_rep=int(open(data_folder+'max_rep').readline().strip())
for number in range(max_rep+1):
	print number
	d=cPickle.load(open(data_folder+'chunk_dict'+str(number), 'rb'))	
	for gene in d:
		if gene not in gene_tracker:
			gene_tracker[gene]=[number]
		else:
			if gene not in special_genes:
				special_genes[gene]={}
			for coord in d[gene]:
				special_genes[gene][coord]=d[gene][coord]
for gene in special_genes:
	d=cPickle.load(open(data_folder+'chunk_dict'+str(gene_tracker[gene][0]), 'rb'))
	for coord in d[gene]:
		special_genes[gene][coord]=d[gene][coord]

spliced_dict={}
for number in range(max_rep+2):
	print number
	if number<max_rep+1:
		d=cPickle.load(open(data_folder+'chunk_dict'+str(number), 'rb'))	
	else:
		d=special_genes
	for gene in d.keys():
		spliced_dict[gene]={}
		tuples=sorted(d[gene].keys())
		for tuple_thing in tuples:
			existing_species=set([])
			length=''
			for species in d[gene][tuple_thing]:
				species_name, species_chrom=species[0].split('.')[:2]
				if not length:
					length=len(species[1])
#				print species_name, good_species[0]
				if species_name==good_species[0]:
					coords=(species_chrom, tuples[0][0], tuples[-1][1])
				existing_species.add(species_name)
			missing_species=species_set-existing_species
			for species in missing_species:
				d[gene][tuple_thing].append([species+'.none', 'N'*length])
			for species_pair in d[gene][tuple_thing]:
				species, sequence=species_pair
				species=species.split('.')[0]
				if coords not in spliced_dict[gene]:
					spliced_dict[gene][coords]={}
				if species not in spliced_dict[gene][coords]:
					spliced_dict[gene][coords][species]=''
				spliced_dict[gene][coords][species]+=sequence.upper()
		if gene in intron_dict:
			for coord in spliced_dict[gene]:
#				print coord
				introns=[]
				for intron in intron_dict[gene]:
					introns.append([intron[0], intron[1]-intron[0]])
				spliced_dict[gene][coord]['intron']=introns
cPickle.dump(spliced_dict, open(data_folder+'spliced_dict', 'wb'), -1)
