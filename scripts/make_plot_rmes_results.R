library(dplyr)
library(ggplot2)
library(ggbeeswarm)
library(stringr)

ptu_info = read.csv('../data/top50ptus_info.csv', header=T, stringsAsFactors = F, row.names = 1)
top50_PTUs = read.csv('../top50-ptus-with-core.txt', header=F, stringsAsFactors = F)$V1

gc_df = read.csv('../output/PTU_core_accessory_GC.csv', header=T, stringsAsFactors = F)

gc_df_diff = data.frame(gc_df %>% group_by(PTU) %>%
  summarise(diff=GC[component=="core"]-GC[component=="accessory"]))
rownames(gc_df_diff) = gc_df_diff$PTU
ptu_info$GC.diff = gc_df_diff[rownames(ptu_info),"diff"]

ptu = "PTU-C"

k = 6

accessory_df = read.csv(paste0('../output/PTU-C_accessory_genes_k', k, '.csv'), header=T, stringsAsFactors = F, row.names = 1)
core_df = read.csv(paste0('../output/PTU-C_core_genes_k', k, '.csv'), header=T, stringsAsFactors = F, row.names = 1)

core_df = core_df[rownames(accessory_df),]


palindromes_k = tolower(read.csv(paste0('../palindromes-k', k, '.txt'), header=F, stringsAsFactors = F)$V1 )

core_df_palindromes = core_df[palindromes_k,]
accessory_df_palindromes = accessory_df[palindromes_k,]

core_df_palindromes$diff.rank.accessory = core_df_palindromes$rank - accessory_df_palindromes$rank
core_df_palindromes$kmer = rownames(core_df_palindromes)
core_df_palindromes$kmer = ordered(core_df_palindromes$kmer, 
                                   levels=core_df_palindromes$kmer[order(core_df_palindromes$diff.rank.accessory)])
# gc of k-mer?
core_df_palindromes$GC.of.kmer = as.factor(sapply(as.character(core_df_palindromes$kmer),
                                        function(x) (stringr::str_count(x, "g")+stringr::str_count(x, "c"))/nchar(x)))
p.avoidance = ggplot(core_df_palindromes, aes(kmer, diff.rank.accessory, fill=GC.of.kmer))+
  geom_bar(stat="identity")+
  scale_fill_brewer(palette="Oranges")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  coord_flip()+
  xlab("")+
  labs(fill="GC content")+
  ylab("Core rank - accessory rank (negative=avoided more in core)")
ggsave(p.avoidance, file='../output/example-PTU-C-k6-avoidance-ranks.pdf')
ggplot(core_df_palindromes, aes(GC.of.kmer, diff.rank.accessory))+
  geom_boxplot()+
  xlab("GC-content of k-mer")+
  ylab("Core-accessory rank (lower=more avoidance in accessory)")

