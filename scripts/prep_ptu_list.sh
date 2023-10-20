while read ptu;
do
	awk -F '\t' '$24=="$ptu"' ../data/table_50PTUs_RS200.csv | awk -F '\t' '{print $1}'> ../data/"$ptu"_accs.txt
	echo $f
done < ../data/PTUs.txt