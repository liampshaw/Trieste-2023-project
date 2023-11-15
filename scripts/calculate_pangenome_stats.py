# From the core/accessory fastas, return pangenome stats
import pandas as pd 
import re

ptu_df_file = 'data/top50ptus_info.csv' # PTU information file (host range)
ptu_info_df = pd.read_csv(ptu_df_file)
dir_prefix = 'core-accessory-fastas/' # prefix for fasta files

def get_headers(input_fasta_file):
	headers = []
	with open(input_fasta_file, 'r') as f:
		for line in f.readlines():
			if line.startswith('>'):
				headers.append(line[1:].strip('\n'))
			else:
				pass
	return(headers)


with open('data/top50ptus_pangenome_stats.csv', 'w') as f:
	f.write('PTU,core_genes_total,core_genes,acc_genes_total,acc_genes\n')
	for ptu in ptu_info_df['PTU']:
		headers = get_headers(dir_prefix+'/'+ptu+'_core_genes.fa')
		total_genes = len(headers)
		gene_groups = len(set([re.sub('.*\\|', '', x) for x in headers]))
		f.write('%s,%d,%d,' % (ptu, total_genes, gene_groups))
		headers = get_headers(dir_prefix+'/'+ptu+'_accessory_genes.fa')
		total_genes = len(headers)
		gene_groups = len(set([re.sub('.*\\|', '', x) for x in headers]))
		f.write('%d,%d\n' % (total_genes, gene_groups))



