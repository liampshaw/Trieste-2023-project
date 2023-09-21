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
