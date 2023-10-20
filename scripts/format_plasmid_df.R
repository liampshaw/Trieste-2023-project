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

write.csv(PTU_df, file='../data/top50ptus_info.csv', row.names = F)
