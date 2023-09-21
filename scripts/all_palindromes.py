
import argparse  
import itertools as iter  
alt_map = {'ins':'0'}
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

def get_options():
    parser = argparse.ArgumentParser()
    parser.add_argument('--N', help='bp length of palindrome', required=True)
    return parser.parse_args()

def reverse_complement(seq):
    for k,v in alt_map.items():
        seq = seq.replace(k,v)
    bases = list(seq)
    bases = reversed([complement.get(base,base) for base in bases])
    bases = ''.join(bases)
    for k,v in alt_map.items():
        bases = bases.replace(v,k)
    return bases

def all_palindromes(N):
    '''Returns all palindromes of length N.'''
    palindromes = []
    for seq in iter.product('ATCG', repeat=N):
        seq = ''.join(seq)
        half_seq = seq[0:int(N/2)]
        end_seq = seq[int(N/2):int(N)]
        if half_seq==reverse_complement(end_seq):
            palindromes.append(seq)
    return palindromes

def main():
    args = get_options()
    palindromes = all_palindromes(int(args.N))
    for p in palindromes:
        print(p)

if __name__=="__main__":
    main()