'''
Run after 'automated_alignments.py'. Strings together scripts to automate
simulation of randomly evolving sequence. Useful if the alignments have been
generated already but more or different simulations are needed. User inputs
data folder where alignments have been generated (same folder where bed file
was input for 'automated_alignments.py') and replicates to generate in a range
from rep_start to rep_end. Simulations will go to a 'simulations' folder.

Dependencies:
1. successfully completed automated_alignments.py or all_together.py

User-configurable parameters:
1. data_folder where bed file is located and output will go
2. species of interest to analyze
3. desired range of replicate numbers to simulate (rep_start inclusive, rep_end
exclusive)
'''


import subprocess
import tree_manips
import time

#Configurable options in this section:
data_folder='/home/alfred/grimson_postdoc/big_data/turnover_rates/example_data/' #tell the program where your bed file is located (output will also go here)
target_species=['hg19', 'nomLeu3', 'calJac3', 'otoGar3', 'canFam3', 'felCat5', 'equCab2', 'oviAri3', 'bosTau7', 'mm10'] #reference genome is mandatory and must be listed first
#common names for above genomes are: 'human', 'gibbon', 'marmoset', 'bushbaby', 'dog', 'cat', 'horse', 'sheep', 'cow', 'mouse'
rep_start, rep_end='0','100'

subprocess.call(['mkdir', stats_folder])
print 'finished preparing utrs, now for simulation'
time.sleep(5)

print 'generating mutation table'
subprocess.call(['python', 'generate_scaled_trinucleotide_table.py', data_folder])

print 'grabbing correct ancestor species number'
reference_ancestor=open('reference_ancestor').readline()

print 'simulating randomly mutated descendants from ancestral sequences, this may take some time'
subprocess.call(['python', 'simulate_di_turnover_trimers.py', data_folder, rep_start, rep_end, target_species[0], reference_ancestor])
