# call like
# ./scripts/count_in_gene.sh
# counts palindrome density at k=6 for each gene group in core/accessory
while read f;
do
	echo $f
	echo "ptu,category,gene_group,n_sequence,average_length,n_palindromes,palindrome_density" > core-accessory-fastas/$f-gene-level-densities.csv
	python scripts/palindrome_density_gene_level.py --fasta core-accessory-fastas/"$f"_core_genes.fa --k 6 --category core --ptu $f >> core-accessory-fastas/$f-gene-level-densities.csv
	 python scripts/palindrome_density_gene_level.py --fasta core-accessory-fastas/"$f"_core_genes.fa --k 6 --category accessory --ptu $f >> core-accessory-fastas/$f-gene-level-densities.csv
done < ptu-list.txt 
