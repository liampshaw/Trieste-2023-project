# Returns total palindrome density of a core/accessory fasta file
# in terms of 

import count_kmers as ck 
import argparse  
from itertools import product
import all_palindromes 
import re

def get_options():
	parser = argparse.ArgumentParser(description='Returns total palindrome density at gene level.',
                                     prog='palindrome_densities')
	parser.add_argument('--fasta', help='fasta file', required=True)
	parser.add_argument('--k', help='value of k', required=True) 
	parser.add_argument('--category', help='category to add to output', required=False,
		default='') 
	parser.add_argument('--ptu', help='PTU name to add to output', required=False,
		default='') 
	return parser.parse_args()


def seq_length(fasta):
	length = 0
	with open(fasta, 'r') as f:
		for line in f.readlines():
			if not line.startswith('>'):
				length += len(line.strip('\n'))
	return(length)


def read_core_or_accessory_fasta(fasta):
	sequences = {}

	current_group = ''
	on_fasta_line = False
	with open(fasta, 'r') as f:
		for line in f.readlines():
			if line.startswith('>'):
				on_fasta_line = False
				group = re.sub('.*\\|', '', line.strip())
				current_group = group
			else:
				if current_group in sequences.keys():
					if on_fasta_line==False:
						sequences[current_group] += 'Z'+line.strip()
						on_fasta_line = True
					elif on_fasta_line==True:
						sequences[current_group] += line.strip()
				else:
					sequences[current_group] = line.strip()
					on_fasta_line=True
	return(sequences)

def main():
	args = get_options()
	input_fasta = args.fasta
	seqs = read_core_or_accessory_fasta(input_fasta)
	#sequence_length = seq_length(input_fasta)
	k = int(args.k)
	for gene_group, sequence in seqs.items():
		sequence_length = len(sequence) - sequence.count('Z')-k + 1# normalised by possible k-mers in whole sequence
		n_sequences_in_group = sequence.count('Z')+1
		output = ck.count_kmers_seq(sequence, k)
		palindromes = all_palindromes.all_palindromes(k)
		kmers = [''.join(kmer) for kmer in product('ACGT', repeat=k)]
		lexicographical_indices_palindromes = [i for i in range(0, 4**k) if kmers[i] in palindromes]
		total_sum = 0
		for x in lexicographical_indices_palindromes:
			total_sum += float(output[x])
		output_string = ','.join([args.ptu, args.category, gene_group, 
			str(n_sequences_in_group),
			str(round(sequence_length/n_sequences_in_group, 1)), str(total_sum), 
			str(round(total_sum/(sequence_length), 6))])
		print(output_string)


if __name__=="__main__":
	main()