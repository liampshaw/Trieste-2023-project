# Check overlap between PTU genera and k-mer database
# previously constructed for previous NAR paper

ptus = read.csv('../data/table_50PTUs_RS200.csv', 
                sep='\t', stringsAsFactors = F)
sort(table(ptus$Genus))

kmer_db = read.csv('~/Downloads/21923121/6-kmer-db.csv',
                   header=T, row.names = 1)
kmer.genera = unique(sapply(colnames(kmer_db),
                     function(x) gsub("_.*", "", x)))
table(kmer.genera %in% names(table(ptus$Genus)))
table(names(table(ptus$Genus)) %in% kmer.genera)

cat(names(table(ptus$Species)))

ptus$Species.simple = gsub("sp..*", "sp.", ptus$Species)

cat(names(table(ptus$Genus)), sep="\n", file="../data/PTU-genera.txt")
