# Compute the density of various targets for plasmids within a PTU

import pandas as pd
import re
# Run from main directory (i.e. call as scripts/compute_target_density_ptu_core_accessory.py)

# Takes:
target_db_file = 'output/target_db.csv'# target database, generated with combined_rmsFinder_output.py
dir_prefix = 'core-accessory-fastas/' # prefix for fasta files
taxonomy_df_file = 'data/PTU-genera-taxonomy.txt' # taxonomy file 
ptu_df_file = 'data/top50ptus_info.csv' # PTU information file (host range)

# need functions for:
# conversion of ambiguous characters YES
# reverse complementing
# reading in a fasta file (and reverse complementing)
# 


alt_map = {'ins':'0'}
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

def reverse_complement(seq):
    for k,v in alt_map.items():
        seq = seq.replace(k,v)
    bases = list(seq)
    bases = reversed([complement.get(base,base) for base in bases])
    bases = ''.join(bases)
    for k,v in alt_map.items():
        bases = bases.replace(v,k)
    return bases

ambiguity_codes =   {'A': ['A'],\
                    'G' : ['G'],\
                    'C'	: ['C'],\
                    'T' : ['T'],\
                    'Y'	: ['C', 'T'],\
                    'R'	: ['A','G'],\
                    'W'	: ['A','T'],\
                    'S'	: ['G','C'],\
                    'K'	: ['T','G'],\
                    'M'	: ['C','A'],\
                    'D'	: ['A','G','T'],\
                    'V'	: ['A','C','G'],\
                    'H'	: ['A','C','T'],\
                    'B'	: ['C','G','T'],\
                    'X'	: ['A','C','G','T'],\
                    'N'	: ['A','C','G','T'],\
                    '-' : ['-']}

def possibleSequences(dna_sequence):
    '''Takes a DNA sequence (possibly containing ambiguous bases) and returns
    a list of all possible sequences.
    Args:
        dna_sequence (str)
            String of DNA sequence (uppercase)
    Returns:
        possible_strings (list)
            List of strings - all possible DNA sequences
    '''
    # Get all possible bases at each position
    possible_bases = [ambiguity_codes[base] for base in dna_sequence]
    # Go through and add these on
    possible_strings = []
    for bases in possible_bases:
        if len(possible_strings)==0: # if no strings yet, use first base
            possible_strings = list(bases)
        elif len(bases)==1: # if just one possible base at position, add it on
            possible_strings = [x+bases[0] for x in possible_strings]
        else: # if multiple possibilities, add them on one by one
            additions = []
            for base in bases:
                additions.append([x+base for x in possible_strings])
            possible_strings = [x for addition in additions for x in addition]
    return possible_strings


def get_seq(input_fasta):
    dna = ''
    with open(input_fasta, 'r') as f:
        for line in f.readlines():
            if line.startswith('>'):
                dna += 'Z'
            else:
                dna += line.strip('\n')
    return(dna)



# ARGUMENTS
# 

# Read in taxonomy df
taxonomy_df = pd.read_csv(taxonomy_df_file, sep=';')
taxonomy_df.drop('Unnamed: 6', axis=1, inplace=True)

# Read in PTU info
ptu_info_df = pd.read_csv(ptu_df_file)
#host_range_map = {'I': 'Species', 'II':'Genus', 'III':'Family', 'IV':'Order', 'V':'Class', 'VI':'Phylum'}
#ptu_info_df['range'] = ptu_info_df['host.range'].map(host_range_map)

# N.B. Treat species as genus (for now) - because e.g. Enterococcus is all Enterococcus faecium in 
# the RefSeq200 top 50 PTUs

ptu="PTU-Y"

for ptu in ptu_info_df['PTU']:
	print(ptu)
	if ptu!='PTU-E9' and ptu!='PTU-E76':
		# Get sequence of core/accessory fastas
		core_seq = get_seq(dir_prefix+ptu+'_core_genes.fa')
		acc_seq = get_seq(dir_prefix+ptu+'_accessory_genes.fa')

		# Target database
		target_df = pd.read_csv(target_db_file)
		target_df['ambiguity'] = [len(possibleSequences(x)) for x in target_df['sequence']]
		target_df['length'] = [len(x) for x in target_df['sequence']]

		# For every unique target, calculate the density in the fasta files


		# Let's just care about non-overlapping cases, so use .count() 
		core_target_densities = {target:sum([core_seq.count(x) + 
									core_seq.count(reverse_complement(x)) 
									for x in possibleSequences(target)])/(2*len(core_seq)-core_seq.count('Z')) 
									for target in set(target_df['sequence'])}
		accessory_target_densities = {target:sum([acc_seq.count(x) + 
									acc_seq.count(reverse_complement(x)) 
									for x in possibleSequences(target)])/(2*len(acc_seq)-acc_seq.count('Z')) 
									for target in set(target_df['sequence'])}

		target_df['core_density'] = target_df['sequence'].map(core_target_densities)
		target_df['acc_density'] = target_df['sequence'].map(accessory_target_densities)

		# Get the host range (included genera) for the PTU
		parent_taxa = ptu_info_df.loc[ptu_info_df['PTU']==ptu]['parent_taxa'].iloc[0]
		host_range = ptu_info_df.loc[ptu_info_df['PTU']==ptu]['range_level'].iloc[0]
		if host_range=='Species':
			parent_taxa = re.sub(' .*', '', parent_taxa)
			host_range = 'Genus'

		genera_within_range = list(taxonomy_df[taxonomy_df[host_range]==parent_taxa]['Genus'])
		# Take the parent taxa and match on that from taxonomy_df using range
		target_df['within_range'] = target_df['genus'].isin(genera_within_range)

		# Write final csv
		target_df.to_csv('output/'+ptu+'-target-counts.csv', index=False)

		# Find all sequences within range
		targets_within_range = set(target_df['sequence'][target_df['within_range']==True])
		targets_without_range = set(target_df['sequence'][target_df['within_range']==False])
		targets_only_within_range = targets_within_range- targets_without_range
		targets_only_without_range = targets_without_range- targets_within_range
		targets_in_both = targets_without_range.intersection(targets_within_range)
		#print('without:', targets_only_without_range)
		#print('within:', targets_only_within_range)
		#print('both:', targets_in_both)

		unique_target_df = target_df.drop_duplicates(subset='sequence', keep='first')
		unique_target_df.drop(['genus', 'genome', 'total_genomes', 'normalized_count'], axis=1, inplace=True)

		def categorize_sequence(seq):
		    if seq in targets_only_within_range:
		        return 'within'
		    elif seq in targets_only_without_range:
		        return 'without'
		    elif seq in targets_in_both:
		        return 'both'
		    else:
		        return 'other'  # Optional: for sequences not in any list

		# Apply the function to the 'sequence' column to create a new column
		unique_target_df['category'] = unique_target_df['sequence'].apply(categorize_sequence)

		#print(unique_target_df)
		unique_target_df.to_csv('output/'+ptu+'-target-counts-unique.csv', index=False)






