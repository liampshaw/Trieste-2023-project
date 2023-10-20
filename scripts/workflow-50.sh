# Make core/accessory gene fastas
while read f;
do
	python scripts/make_core_accessory_fasta_from_roary.py \
		--ffn_dir /Users/Liam/Downloads/ffn \
		--gene_presence_absence /Users/Liam/Downloads/roary/"$f"_roary_i80_c80/gene_presence_absence.Rtab \
		--clustered_proteins /Users/Liam/Downloads/roary/"$f"_roary_i80_c80/clustered_proteins \
		--core 0.8 \
		--output core-accessory-fastas/"$f"
	echo $f 
done < top50-ptus.txt 

# Calculate palindrome densities
while read f;
do
	# Repeat palindrome density analysis for these fastas
	python scripts/palindrome_densities.py --k 4 --fasta core-accessory-fastas/"$f"_core_genes.fa > core-accessory-palindromes/"$f"_k4_core_palindrome_densities.txt
	python scripts/palindrome_densities.py --k 4 --fasta core-accessory-fastas/"$f"_accessory_genes.fa > core-accessory-palindromes/"$f"_k4_palindrome_densities.txt
	python scripts/palindrome_densities.py --k 6 --fasta core-accessory-fastas/"$f"_core_genes.fa > core-accessory-palindromes/"$f"_k6_core_palindrome_densities.txt
	python scripts/palindrome_densities.py --k 6 --fasta core-accessory-fastas/"$f"_accessory_genes.fa > core-accessory-palindromes/"$f"_k6_accessory_palindrome_densities.txt
	echo $f
done < top50-ptus.txt 


# GC contents
while read ptu;
do
	echo $ptu,"accessory",$(python scripts/GC_content.py --fasta core-accessory-fastas/"$ptu"_accessory_genes.fa)
	echo $ptu,"core",$(python scripts/GC_content.py --fasta core-accessory-fastas/"$ptu"_core_genes.fa)
done < top50-ptus.txt > output/PTU_core_accessory_GC.csv
