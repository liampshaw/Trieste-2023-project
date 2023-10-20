#!/bin/bash
# for running on BMRC cluster (Liam Shaw)
for f in *fa;
do
	name=$(echo $f | cut -d '.' -f 1)
	sbatch -p short run-rmes.sh $f "$name"_k4.csv 4
	sbatch -p short run-rmes.sh $f "$name"_k6.csv 6
done 
