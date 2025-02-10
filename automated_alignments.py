'''
Strings together several scripts to automatically generate alignments 
(including reconstructing ancestral sequences). User inputs data folder where
an input bed file is located, where output data will also go, a folder where
maf files (organized by chromosome) can be found, a desired subset of species
names to be examined (reference species first), the name of the input bed file,
and the name of an input phylogeny topology to use.

Dependencies:
1. ete2 python package in a unix $PYTHONPATH folder (from etetoolkit.org)
2. dnaml binary in a unix $PATH folder (from phylip)
3. downloaded maf files sorted by a reference species (mine are from UCSC)
4. newick file showing how species in maf files are related
5. bed file from same reference genome as maf file with UCSC formatting and one
record per gene, carefully filtered to avoid redundancies and bad genes

User-configurable parameters:
1. data_folder where bed file is located and output will go
2. name of bed file
3. folder where maf files are stored
4. species of interest to analyze
5. name of newick formatted file (must be stored in current folder)
6. name of output folder where stats will go (created within the data folder)
'''

import subprocess
import tree_manips
import time

#Configurable options in this section:
data_folder='/home/alfred/grimson_postdoc/big_data/turnover_rates/example_data/' #tell the program where your bed file is located (output will also go here)
input_bed_file='4-28-15_refseq.bed' #replace this with your own bed file
maf_folder='/media/alfred/data_drive/big_data/reference_files/multiz_100way/' #tell the program where maf files are located
target_species=['hg19', 'nomLeu3', 'calJac3', 'otoGar3', 'canFam3', 'felCat5', 'equCab2', 'oviAri3', 'bosTau7', 'mm10'] #reference genome is mandatory and must be listed first
#common names for above genomes are: 'human', 'gibbon', 'marmoset', 'bushbaby', 'dog', 'cat', 'horse', 'sheep', 'cow', 'mouse'
input_tree='hg19.100way.nh' #replace this with your own newick file
stats_folder=data_folder+'output_stats_folder/' #you can change output_stats_folder to something more memorable if desired

#main programs in this section:
subprocess.call(['mkdir', stats_folder])
bed_prefix=input_bed_file.replace('.bed', '')

print 'converting...'
subprocess.call(['python', 'coordinate_converter.py', data_folder, input_bed_file])

print 'getting longest isoform...'
subprocess.call(['python', 'get_longest_UTRs.py', data_folder, bed_prefix+'_3UTR_exons.bed'])

print 'splitting and sorting bed files by chromosome...'
subprocess.call(['python', 'prepare_utrs.py', data_folder, bed_prefix+'_3UTR_exons_best_isoform.bed'])

print 'extracting maf overlaps'
subprocess.call(['python', 'find_maf_overlaps.py', data_folder, maf_folder, str(target_species)])

print 'adding missing species, splicing introns...'
subprocess.call(['python', 'add_species_rm_introns.py', data_folder, str(target_species), bed_prefix])

print 'removing ambiguous nucleotides, reverse complementing'
subprocess.call(['python', 'remove_gaps.py', data_folder, str(target_species), str(target_species), bed_prefix+'_3UTR_exons_best_isoform', 'spliced_dict'])

print 'naming internal nodes of phylogenetic tree...'
subprocess.call(['python', 'name_internal_nodes.py', data_folder, input_tree])

print 'pruning tree and renaming for dnaml ancestor reconstruction...'
tree_manips.prune_species(data_folder+input_tree.replace('.nh', '')+'mod.nh', target_species, 'intree')

print 'reconstructing ancestors...'
subprocess.call(['python', 'make_ancestors.py', data_folder, str(target_species)])

print 'reformatting, uppercasing...'
subprocess.call(['python', 'reformat_utr_dict.py', data_folder, str(target_species)])
