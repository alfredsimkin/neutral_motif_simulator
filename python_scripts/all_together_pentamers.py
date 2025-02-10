'''
Unlike all_together, this is set up to model pentamers: instead of looking at
trinucleotides whose middle nucleotide is free to change, looks at
pentanucleotides whose middle nucleotide is free to change.

Strings together several scripts to automate analysis from start to finish.
User inputs data folder where an input bed file is located, where output data
will also go, a folder where maf files (organized by chromosome) can be found,
a desired subset of species names to be examined (reference species first), the
name of the input bed file, the name of an input phylogeny topology to use, and
replicates to analyze in a range from rep_start to rep_end. Simulations will go
to a 'simulations' folder, and final output will go to a 'stats' folder, both
within the data folder

Dependencies:
1. ete2 python package in a unix $PYTHONPATH folder (from etetoolkit.org)
2. matplotlib python package
3. scipy python package
4. dnaml binary in a unix $PATH folder (from phylip)
5. downloaded and unzipped maf files sorted by a reference species (mine are from UCSC)
6. newick file showing how species in maf files are related
7. bed file from same reference genome as maf file with UCSC formatting and one
record per gene, carefully filtered to avoid redundancies and bad genes
8. file containing miRNAs (or other 8mers) of interest for statistical analysis
(generated by populate_real_mers.py)

User-configurable parameters:
1. data_folder where bed file is located and output will go
2. name of bed file
3. folder where maf files are stored
4. species of interest to analyze
5. name of newick formatted file (must be stored in current folder)
6. name of output folder where stats will go (created within the data folder)
7. full path to file with miRNAs to be analyzed
8. desired range of replicate numbers to analyze (rep_start inclusive, rep_end
exclusive)
'''


import subprocess
import tree_manips
import time

#Configurable options in this section:
data_folder='/mnt/hulk_tank/alfred_stuff/silver_surfer/Grimson_new_phylogeny/data/pentamers/' #tell the program where your bed file is located (output will also go here)
input_bed_file='some_bed_file' #replace this with your own bed file
maf_folder='/mnt/hulk_tank/alfred_stuff/silver_surfer/Grimson_new_phylogeny/data/maf_data_hg38_100_way/' #tell the program where maf files are located
target_species=['hg38', 'oryCun2', 'nomLeu3', 'susScr3', 'calJac3', 'canFam3', 'macFas5', 'bosTau8', 'mm10', 'monDom5'] #reference genome is mandatory and must be listed first
#common names for above genomes are: Human, Rabbit, Gibbon, Pig, Marmoset, Dog, Crab-eating macaque, Cow, Mouse, Opossum
input_tree='hg38.100way.nh' #replace this with your own newick file
stats_folder=data_folder+'stats/' #you can change output_stats_folder to something more memorable if desired
miRNA_file='/mnt/hulk_tank/alfred_stuff/silver_surfer/Grimson_new_phylogeny/data/real_eightmers_human_mouse_fish_best'
rep_start, rep_end='0','100'

#main body in this section:
subprocess.call(['mkdir', stats_folder])
bed_prefix=input_bed_file.replace('.bed', '')

print 'converting...'
#subprocess.call(['python', 'coordinate_converter.py', data_folder, input_bed_file])

print 'getting longest isoform...'
#subprocess.call(['python', 'get_longest_UTRs.py', data_folder, bed_prefix+'_3UTR_exons.bed'])

print 'splitting and sorting bed files by chromosome...'
#subprocess.call(['python', 'prepare_utrs.py', data_folder, bed_prefix+'_3UTR_exons_best_isoform.bed'])

print 'extracting maf overlaps'
#subprocess.call(['python', 'find_maf_overlaps.py', data_folder, maf_folder, str(target_species)])

print 'adding missing species, splicing introns...'
#subprocess.call(['python', 'add_species_rm_introns.py', data_folder, str(target_species), bed_prefix])

print 'removing ambiguous nucleotides, reverse complementing'
#subprocess.call(['python', 'remove_gaps.py', data_folder, str(target_species), str(target_species), bed_prefix+'_3UTR_exons_best_isoform', 'spliced_dict'])

print 'naming internal nodes of phylogenetic tree...'
#subprocess.call(['python', 'name_internal_nodes.py', data_folder, input_tree])

print 'pruning tree and renaming for dnaml ancestor reconstruction...'
#tree_manips.prune_species(data_folder+input_tree.replace('.nh', '')+'mod.nh', target_species, 'intree')

print 'reconstructing ancestors...'
#subprocess.call(['python', 'make_ancestors.py', data_folder, str(target_species)])

print 'reformatting, uppercasing...'
#subprocess.call(['python', 'reformat_utr_dict.py', data_folder, str(target_species)])

print 'finished preparing utrs, now for simulation'
time.sleep(5)

print 'generating mutation table'
subprocess.call(['python', 'generate_scaled_pentanucleotide_table.py', data_folder])

print 'grabbing correct ancestor species number'
reference_ancestor=open('reference_ancestor').readline()

print 'simulating randomly mutated descendants from ancestral sequences, this may take some time'
subprocess.call(['python', 'simulate_di_turnover_pentamers.py', data_folder, rep_start, rep_end, target_species[0], reference_ancestor])

print 'counting real turnover events'
subprocess.call(['python', 'count_real.py', data_folder, stats_folder])

print 'making distributions'
subprocess.call(['python', 'make_distributions_expanded.py', stats_folder, rep_start, rep_end])

#ranks miRNA gain and loss rates relative to other 8mers and send to files
print 'ranking miRNAs'
subprocess.call(['python', 'extract_miRNA_ranks.py', stats_folder, miRNA_file])

#plots the rates from extract_miRNA_ranks.py
print 'plotting miRNA ranks'
subprocess.call(['python', 'gains_vs_losses.py', stats_folder])

#bins all 8mers and conserved miRNAs by their observed number of std devs away
#from simulated values, calculates stats for each bin
print 'binning all 8mers and miRNA 8mers into stdev bins'
subprocess.call(['python', 'bin_distributions.py', stats_folder])

print 'graphing the bins'
subprocess.call(['python', 'graph_distributions_matplotlib.py', stats_folder])
