import custom
import cPickle
import subprocess
import sys
data_folder, species_list=sys.argv[1:]
species_list=eval(species_list)
reference=species_list[0]
total_species=len(species_list)
utr_dict=cPickle.load(open(data_folder+'spliced_dict_nogaps', 'rb'))
ancestor_file=open('reference_ancestor', 'w')
dnaml_output=open('dnaml_input', 'w')
dummy_outfile=open('outfile', 'w')
dummy_outtree=open('outtree', 'w')
dummy_infile=open('infile', 'w')
dummy_outfile.close()
dummy_outtree.close()
dummy_infile.close()

def gather_ancestors(infile):
	new_dict, first_line={}, 'big'
	for line_number, line in enumerate(open(infile)):
		line=line.strip()
		if line.startswith('Between'):
			first_line=line_number+3
		if line_number>=first_line:
			line=line.split()
			if len(line)>1:
				new_dict[line[1]]=line[0]
			else:
				break
	return new_dict



phylip_species_list=[]
for species_number, species in enumerate(species_list):
	phylip_species_list.append(species+' '*(10-len(species)))

dnaml_output.write('r\nu\n5\no\n'+str(total_species)+'\ny\nr')
dnaml_output.close()

ancestor_dict={}
fasta_list=[]
not_gathered=True
for gene in utr_dict:
	alignment_list=[]
	for species in species_list:
		alignment_list.append(utr_dict[gene][species])
		fasta_list.append([gene+'_'+species, utr_dict[gene][species]])
	custom.print_phylip(phylip_species_list, [gene, alignment_list], 'infile')
	subprocess.call(['sh', 'dnaml_script.sh'])
	first_line='big'
	if not_gathered==True:
		relation_dict=gather_ancestors('outfile')
		reference_ancestor=relation_dict[reference]
		print relation_dict
		cPickle.dump(relation_dict, open(data_folder+'ancestor_dict', 'w'))
		not_gathered=False
	for line_number, line in enumerate(open('outfile')):
		if line.startswith('Probable'):
			first_line=line_number+4
		if line_number>=first_line:
			line=line.split()
			if len(line)>1:
				species, sequence=line[0], ''.join(line[1:])
				if gene not in ancestor_dict:
					ancestor_dict[gene]={}
				if species not in ancestor_dict[gene]:
					ancestor_dict[gene][species]=''
				ancestor_dict[gene][species]+=sequence
#custom.print_fasta(fasta_list, data_folder+'fasta_formatted_genes', 'w')
ancestor_file.write(reference_ancestor)
cPickle.dump(ancestor_dict, open(data_folder+'anc_aln', 'w'))
cPickle.dump(relation_dict, open(data_folder+'ancestor_dict', 'w'))
