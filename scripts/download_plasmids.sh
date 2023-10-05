for file in *accs.txt;
do
	echo $file
	ptu=$(echo $file | cut -d '_' -f 1)
	mkdir -p $ptu
	while read acc;
	do
		ncbi-acc-download -F fasta $acc
		mv "$acc".fa $ptu/
	done < $file
done 