# Returns GC content of a fasta file

import count_kmers as ck 
import argparse  
from itertools import product
import all_palindromes 

def get_options():
	parser = argparse.ArgumentParser(description='Returns GC content of a fasta file.',
                                     prog='GC_content')
	parser.add_argument('--fasta', help='fasta file', required=True)
	return parser.parse_args()


def GC(fasta):
	seq = ''
	with open(fasta, 'r') as f:
		for line in f.readlines():
			if not line.startswith('>'):
				seq += line.strip('\n')
	gc_content = (seq.count('G')+seq.count('g')+seq.count('C')+seq.count('c'))/len(seq)
	return(gc_content)

if __name__=="__main__":
	args = get_options()
	print(GC(args.fasta))