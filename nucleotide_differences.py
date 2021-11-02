#!/usr/bin/env python3


#run using: python nucleotide_differences.py input.fasta

from itertools import combinations
import sys
import pandas as pd

with open(sys.argv[1]) as f_input, open('output_oneline.fasta', 'w') as f_output:
    block = []

    for line in f_input:
        if line.startswith('>'):
            if block:
                f_output.write(''.join(block) + '\n')
                block = []
            f_output.write(line)
        else:
            block.append(line.strip())

    if block:
        f_output.write(''.join(block) + '\n')

fasta = {}
with open('output_oneline.fasta',) as file_one:
    for line in file_one:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            active_sequence_name = line[1:]
            if active_sequence_name not in fasta:
                fasta[active_sequence_name] = []
            continue
        sequence = line
        fasta[active_sequence_name].append(sequence)

dna_dict = {k:list(v[0]) for k,v in fasta.items()}
substitutions_dict = {}

def PairwiseComparison(d1, d2, dna1, dna2):
	#convert each set of lists into a pandas dataframe
	d1 = pd.DataFrame(d1)
	d2 = pd.DataFrame(d2)
	
	#get number of n's and -'s in the sequences
	n_d1_filter = d1.loc[d1[0] == 'n']
	n_d1 = len(n_d1_filter)
	n_d2_filter = d2.loc[d2[0] == 'n']
	n_d2 = len(n_d2_filter)
	dash_d1_filter = d1.loc[d1[0] == '-']
	dash_d1 = len(dash_d1_filter)
	dash_d2_filter = d2.loc[d2[0] == '-']
	dash_d2 = len(dash_d2_filter)
	
	#reset the index so we have a record of which position each 
	#nucleotide is in
	d1.reset_index(level=0, inplace=True)
	d2.reset_index(level=0, inplace=True)
	
	#filter out missing data from both sequences
	d1 = d1.rename(columns = {d1.columns[1]: 'alleles'})
	d2 = d2.rename(columns = {d2.columns[1]: 'alleles'})
	
	d1_filter = d1.drop(d1[d1.alleles.str.contains(r'[n]') | d1.alleles.str.contains(r'[-]')].index)
	d2_filter = d2.drop(d2[d2.alleles.str.contains(r'[n]') | d2.alleles.str.contains(r'[-]')].index)
	
	#merge the filtered sequences based on index -- this preserves indices from fasta file
	#anything that is different due to missing data is now removed from the dataframe
	paired = pd.merge(d1_filter, d2_filter, on = 'index')
	
	#split them back up to compare (probably a more efficient way to do this but I'm feeling lazy)
	d1_new = pd.DataFrame(paired['alleles_x'])
	d2_new = pd.DataFrame(paired['alleles_y'])
	
	#rename the columns so they'll compare nicely
	d1_new = d1_new.rename(columns = {d1_new.columns[0]: 'alleles'})
	d2_new = d2_new.rename(columns = {d2_new.columns[0]: 'alleles'})
	
	#now compare the two sequences and count the number of subsitutes
	subs = d2_new[d1_new.ne(d2_new).any(axis=1)]
	#get number of substitutions
	substitutions = len(subs)

	print('Number of n in first sequence: ' + str(n_d1), flush = True)
	print('Number of n in second sequence: ' + str(n_d2), flush = True)
	print('Number of - in first sequence: ' + str(dash_d1), flush = True)
	print('Number of - in second sequence: ' + str(dash_d1), flush = True)
	print('Total nucleotide substitutions: ' + str(substitutions), flush = True)
	
	#import substitutions dict for outputting matrix
	global substitutions_dict
	
	pair = str(dna1) + ':' + str(dna2)
	substitutions_dict[pair] = str(substitutions)
	
	
# Creates all the possible combinations of the dictionary keys
dna_combinations = combinations(dna_dict, 2)
sys.stdout = open('compare.log', 'w')

for dna1, dna2 in dna_combinations:
	print(f'first sequence: {dna1} vs second sequence: {dna2}', flush = True)
	comparison = PairwiseComparison(dna_dict[dna1], dna_dict[dna2], dna1, dna2)

sys.stdout = sys.__stdout__

def output_matrix():
	#read back in the filled dictionary with names of population pairs and substitutions
	global substitutions_dict
	
	#convert to a dataframe
	subs = pd.DataFrame(substitutions_dict.items(), columns = ['pairs','substitutions'])
	#split the population pairs back into different columns
	subs['pop1'] = subs['pairs'].str.split(':').map(lambda x: x[0])
	subs['pop2'] = subs['pairs'].str.split(':').map(lambda x: x[1])
	#convert substitutions to numeric
	subs['substitutions'] = pd.to_numeric(subs['substitutions'])
	subs = subs.pivot(index = 'pop1', columns = 'pop2', values = 'substitutions')
	subs.index.name = None
	subs.columns.name = None
	subs.to_csv('substitution_matrix.csv', sep = ',')

output_matrix()
	
	
	
	
