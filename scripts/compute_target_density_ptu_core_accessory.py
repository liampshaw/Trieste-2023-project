# Compute the density of various targets for plasmids within a PTU

# Run from main directory (i.e. call as scripts/compute_target_density_ptu_core_accessory.py)

# Takes:
# target database (output/target_db.csv) done
# core/accessory fasta prefix for PTU
ptu = 'PTU-Y' # name of PTU
dir_prefix = 'core-accessory-fastas/PTU-Y_' # prefix
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
host_range_map = {'I': 'Species', 'II':'Genus', 'III':'Family', 'IV':'Order', 'V':'Class', 'VI':'Phylum'}
ptu_info_df['range'] = ptu_info_df['host.range'].map(host_range_map)

# Treat species as genus (for now) - because e.g. Enterococcus is all Enterococcus faecium in 
# the RefSeq200 top 50 PTUs




# Get sequence of core/accessory fastas
core_seq = get_seq(dir_prefix+ptu+'_core_genes.fa')
core_seq_rev_comp = reverse_complement(core_seq)
acc_seq = get_seq(dir_prefix+ptu+'_accessory_genes.fa')
acc_seq_rev_comp = reverse_complement(core_seq)


# Let's just care about non-overlapping cases, so use .count() 
occurrences = sum([core_seq.count(x) for x in possibleSequences(target)])





