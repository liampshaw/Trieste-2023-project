# Convert pangenome files into core and accessory gene files

import argparse  
import pandas as pd
from math import floor

def get_options():
	parser = argparse.ArgumentParser(description='Creates core and accessory gene files.',
                                     prog='make_core_accessory_fasta')
	parser.add_argument('--pan_genome_reference', help='pan_genome_reference.fa fasta file (roary output)', required=True)
	parser.add_argument('--gene_presence_absence', help='pan_genome_reference.Rtab (roary output)', required=True) 
	parser.add_argument('--core', help='threshold for core (default=0.8)', required=False, default=0.8) 
	parser.add_argument('--output', help='output prefix', required=False, default='output') 
	return parser.parse_args()

def read_pan_genome_reference_fasta(fasta):
	seq_dict = {}
	current_family = ''
	for line in open(fasta, 'r').readlines():
		if line.startswith('>'):
			family_name = line.strip('\n').split(' ')[1]
			seq_dict[family_name] = ''
			current_family = family_name
		else:
			seq = line.strip('\n')
			seq_dict[current_family] += seq
	return(seq_dict)


def main():
	args = get_options()
	gene_sequences = read_pan_genome_reference_fasta(args.pan_genome_reference)
	df = pd.read_csv(args.gene_presence_absence, sep='\t', index_col=0)
	number_of_genomes = len(df.columns)
	core_threshold_n = floor(number_of_genomes*args.core)
	print(core_threshold_n)

	row_sums = df.sum(axis=1).to_dict()
	core_genes = [x for x in row_sums.keys() if row_sums[x]>core_threshold_n]
	accessory_genes = [x for x in row_sums.keys() if row_sums[x]<=core_threshold_n]

	with open(args.output+'_core_genes.fa', 'w') as f:
		for gene in core_genes:
			f.write('>%s\n%s\n' % (gene, gene_sequences[gene]))
	with open(args.output+'_accessory_genes.fa', 'w') as f:
		for gene in accessory_genes:
			f.write('>%s\n%s\n' % (gene, gene_sequences[gene]))


if __name__=="__main__":
	main()