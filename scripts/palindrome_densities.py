# Returns palindrome density of a fasta file in table format

import count_kmers as ck 
import argparse  
from itertools import product
import all_palindromes 

def get_options():
	parser = argparse.ArgumentParser(description='Returns palindrome density of a fasta file.',
                                     prog='palindrome_densities')
	parser.add_argument('--fasta', help='fasta file', required=True)
	parser.add_argument('--k', help='value of k', required=True) 
	return parser.parse_args()


def seq_length(fasta):
	length = 0
	with open(fasta, 'r') as f:
		for line in f.readlines():
			if not line.startswith('>'):
				length += len(line.strip('\n'))
	return(length)


def main():
	args = get_options()
	input_fasta = args.fasta
	sequence_length = seq_length(input_fasta)
	k = int(args.k)
	output = ck.count_kmers(input_fasta, k)
	palindromes = all_palindromes.all_palindromes(k)

	kmers = [''.join(kmer) for kmer in product('ACGT', repeat=k)]
	lexicographical_indices_palindromes = [i for i in range(0, 4**k) if kmers[i] in palindromes]

	for x in lexicographical_indices_palindromes:
		print(kmers[x], float(output[x])/sequence_length)



if __name__=="__main__":
	main()