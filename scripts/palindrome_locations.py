import all_palindromes 
import argparse  
import count_kmers


def get_options():
	parser = argparse.ArgumentParser(description='Returns palindrome locations in a fasta file.',
                                     prog='palindrome_locations')
	parser.add_argument('--fasta', help='fasta file', required=True)
	parser.add_argument('--k', help='value of k', required=True) 
	return parser.parse_args()


def palindrome_locations(sequence, k):
	palindromes = all_palindromes.all_palindromes(k)
	# loop through sequence with sliding window of k and count if palindrome or not
	locations = []
	for i in range(0, len(sequence)-k):
		kmer = sequence[i:(i+k)]
		if kmer in palindromes:
			locations.append(i)
	return(locations)


def main():
	args = get_options()
	input_fasta = args.fasta
	seq = count_kmers.get_seq(input_fasta)[1:] # we strip off a Z that is added by default
	k = int(args.k)

	palindromes = all_palindromes.all_palindromes(k)

	locations = palindrome_locations(seq, k)
	for l in locations:
		print(l)

if __name__=="__main__":
	main()


