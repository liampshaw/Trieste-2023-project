# try to format plasmid dataframe to get PTU info
library(dplyr)
plasmid_df = read.csv('../data/table_50PTUs_RS200.csv', header=T, stringsAsFactors = F,
                      sep='\t')

plasmid_df = plasmid_df[order(plasmid_df$PTU),]

# Summarise at level of PTU
PTU_df = plasmid_df %>% group_by(PTU) %>%
  summarise(host.range=unique(HRange),
            median.length=median(Size))

PTU_df$PTU[PTU_df$PTU=="PTU-B/O/K/Z"] = "PTU-B_O_K_Z"
PTU_df$PTU[PTU_df$PTU=="PTU-L/M"] = "PTU-L_M"
plasmid_df$PTU[plasmid_df$PTU=="PTU-B/O/K/Z"] = "PTU-B_O_K_Z"
plasmid_df$PTU[plasmid_df$PTU=="PTU-L/M"] = "PTU-L_M"

host_range_map = c("Species", "Genus", "Family", "Order", "Class", "Phylum")
names(host_range_map) = c("I", "II", "III", "IV", "V","VI")

PTU_df$range_level = as.character(sapply(PTU_df$host.range, 
                            function(x) host_range_map[x]))
PTU_df$parent_taxa = sapply(1:nrow(PTU_df),
                            function(x) head(plasmid_df[which(plasmid_df$PTU==PTU_df$PTU[x]),PTU_df$range_level[x]], 1))

write.csv(PTU_df, file='../data/top50ptus_info.csv', row.names = F)
