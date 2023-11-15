library(ggplot2)
library(tidyr)
library(ggsignif)
library(dplyr)

ptu_db = read.csv('data/top50ptus_info.csv', header=T)

ptu_pangenome_stats = read.csv('data/top50ptus_pangenome_stats.csv', header=T)
ptu_db = merge(ptu_db, ptu_pangenome_stats, by="PTU")
ptu_list = ptu_db$PTU
ptu_list = ptu_list[-which(ptu_list=="PTU-E76")]
ptu_list = ptu_list[-which(ptu_list=="PTU-E9")]

for (ptu in ptu_list){
  df = read.csv(paste0('output/', ptu, '-target-counts-unique.csv'), header=T)

  
  long_df <- pivot_longer(df, 
                          cols = c(acc_density, core_density),
                          names_to = "variable",
                          values_to = "density")
  # deambiguate (normalise by ambiguity)
  long_df$density_deambig = long_df$density/long_df$ambiguity
  
  # Subset to k=6
  long_df_6 = long_df[which(long_df$length==6),]
  
  # Summarise the number of targets
  counts_for_categories = c(paste0("both\n(n=", length(which(long_df_6$category=="both")), ")"),
                            paste0("within\n(n=",  length(which(long_df_6$category=="within")), ")"),
                            paste0("without\n(n=",  length(which(long_df_6$category=="without")), ")"))
  names(counts_for_categories) = c("both", "within", "without")
  long_df_6$category = counts_for_categories[long_df_6$category]
  
  p = ggplot(long_df_6, aes(category, density_deambig, colour=variable))+
    geom_boxplot()+
    xlab("observed RM systems\nin relation to host range of PTU")+
    ylab("density (normalized by ambiguity)")+
    theme_bw()+
    theme(panel.grid = element_blank())+
    ggtitle(paste0(ptu, " (", as.character(round(ptu_db[which(ptu_db$PTU==ptu), "median.length"]/1000, 1)), "kb)", 
                   "\n", "Host range: ", 
                   ptu_db[which(ptu_db$PTU==ptu), "host.range"],
                   " (", ptu_db[which(ptu_db$PTU==ptu), "parent_taxa"], ")\n",
                   as.character(ptu_db[which(ptu_db$PTU==ptu), "core_genes"]), " hard shell genes (n=", 
                   as.character(ptu_db[which(ptu_db$PTU==ptu), "core_genes_total"]), " total); ",
                   as.character(ptu_db[which(ptu_db$PTU==ptu), "acc_genes"]), " accessory genes (n=", 
                   as.character(ptu_db[which(ptu_db$PTU==ptu), "acc_genes_total"]), " total)"))
  ggsave(file=paste0('output/target-plots/', ptu, '_k-6.pdf'),
         p)
  
}


#   ggsignif::geom_signif(comparisons=list(c("acc_density", "core_density")),
#                           method="wilcox")
# 
# # df$diff = df$core_density-df$acc_density
# 
# # long_df_6_mutate = long_df_6 %>% group_by(sequence) %>%
#   mutate(diff=density_deambig[variable=="core_density"]-density_deambig[variable=="acc_density"])
# ggplot(long_df_6_mutate, aes(diff, fill=category))+
#   geom_histogram()+
#   facet_wrap(~category, nrow=3)


