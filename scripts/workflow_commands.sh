#!/bin/bash

# Palindrome locations (for whole plasmid)
for f in PTU-*/*.fa;
do
	PTU=$(echo $f | cut -d '/' -f 1)
	plasmid=$(echo $f | cut -d '/' -f 2 | rev | cut -d '.' -f 2- | rev)
	python scripts/palindrome_locations.py --k 4 --fasta $f > "$PTU"/"$plasmid"_k4_palindrome_locations.txt
	python scripts/palindrome_locations.py --k 6 --fasta $f > "$PTU"/"$plasmid"_k6_palindrome_locations.txt
	echo $f
done

# Palindrome densities (should also work on multi-fasta files i.e. core genes vs. accessory genes.)
for f in PTU-*/*.fa;
do
	PTU=$(echo $f | cut -d '/' -f 1)
	plasmid=$(echo $f | cut -d '/' -f 2 | rev | cut -d '.' -f 2- | rev)
	python scripts/palindrome_densities.py --k 4 --fasta $f > "$PTU"/"$plasmid"_k4_palindrome_densities.txt
	python scripts/palindrome_densities.py --k 6 --fasta $f > "$PTU"/"$plasmid"_k6_palindrome_densities.txt
	echo $f
done

# Make core/accessory gene fastas
for ptu in PTU-*/;
do
	python scripts/make_core_accessory_fasta.py --pan_genome_reference "$ptu"/pan_genome_reference.fa --gene_presence_absence "$ptu"/gene_presence_absence.Rtab \
		--output "$ptu"/component
	# Repeat palindrome density analysis for these fastas
	python scripts/palindrome_densities.py --k 4 --fasta "$ptu"/component_core_genes.fa > "$ptu"/k4_core_palindrome_densities.txt
	python scripts/palindrome_densities.py --k 4 --fasta "$ptu"/component_accessory_genes.fa > "$ptu"/k4_accessory_palindrome_densities.txt
	python scripts/palindrome_densities.py --k 6 --fasta "$ptu"/component_core_genes.fa > "$ptu"/k6_core_palindrome_densities.txt
	python scripts/palindrome_densities.py --k 6 --fasta "$ptu"/component_accessory_genes.fa > "$ptu"/k6_accessory_palindrome_densities.txt
done

