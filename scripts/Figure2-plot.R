# Figure 2
# Design: example broad host-range plasmid with normalized density

library(ggplot2)
library(tidyr)
library(ggsignif)
library(dplyr)

ptu_db = read.csv('../data/top50ptus_info.csv', header=T)

ptu_pangenome_stats = read.csv('../data/top50ptus_pangenome_stats.csv', header=T)
ptu_db = merge(ptu_db, ptu_pangenome_stats, by="PTU")
ptu_list = ptu_db$PTU
ptu_list = ptu_list[-which(ptu_list=="PTU-E76")]
ptu_list = ptu_list[-which(ptu_list=="PTU-E9")]

for (ptu in ptu_list){
  df = read.csv(paste0('../output/', ptu, '-target-counts-unique.csv'), header=T)
  
  
  long_df <- pivot_longer(df, 
                          cols = c(acc_density, core_density),
                          names_to = "variable",
                          values_to = "density")
  # deambiguate (normalise by ambiguity)
  long_df$density_deambig = long_df$density/long_df$ambiguity
  
  # Subset to k=6
  long_df_6 = long_df[which(long_df$length==6),]
  
  # Summarise the number of targets
  counts_for_categories = c(paste0("ubiquitous targets\n(n=", length(unique(long_df_6$sequence[which(long_df_6$category=="both")])), ")"),
                            paste0("targets inside host range\n(n=",  length(unique(long_df_6$sequence[which(long_df_6$category=="within")])), ")"),
                            paste0("targets outside host range\n(n=",  length(unique(long_df_6$sequence[which(long_df_6$category=="without")])), ")"))
  names(counts_for_categories) = c("both", "within", "without")
  long_df_6$category = counts_for_categories[long_df_6$category]
  
  long_df_6$variable = sapply(long_df_6$variable,
                              function(x) ifelse(x=="acc_density", "accessory", "hard shell"))
  p = ggplot(long_df_6, aes(density_deambig, x=variable, colour=variable))+
    geom_hline(yintercept = 1/4096, linetype='dashed')+
    geom_boxplot()+
    facet_wrap(~category, nrow=1)+
    xlab("")+
    ylab("density (normalized by ambiguity)")+
    theme_bw()+
    scale_color_manual(values=c("black", "red"))+
    theme(panel.grid = element_blank())+
    ggtitle(paste0(ptu, " (", as.character(round(ptu_db[which(ptu_db$PTU==ptu), "median.length"]/1000, 1)), "kb)", 
                   "\n", "Host range: ", 
                   ptu_db[which(ptu_db$PTU==ptu), "host.range"],
                   " (", ptu_db[which(ptu_db$PTU==ptu), "parent_taxa"], ")\n",
                   as.character(ptu_db[which(ptu_db$PTU==ptu), "core_genes"]), " hard shell genes (n=", 
                   as.character(ptu_db[which(ptu_db$PTU==ptu), "core_genes_total"]), " total)\n",
                   as.character(ptu_db[which(ptu_db$PTU==ptu), "acc_genes"]), " accessory genes (n=", 
                   as.character(ptu_db[which(ptu_db$PTU==ptu), "acc_genes_total"]), " total)"))+
    theme(legend.position = "none")+
    theme(axis.text=element_text(colour="black"))+
    ggsignif::geom_signif(comparisons=list(c("accessory", "hard shell")),
                          test = "wilcox.test")
  assign(paste0("p.", gsub("PTU-", "", ptu)), p) 
  ggsave(p, file=paste0("../figures/", Sys.Date(), "-", ptu, ".pdf"))
}


pdf(paste0('../figures/', Sys.Date(), '-Figure-2.pdf'), width=14, height=12)
p.small = cowplot::plot_grid(p.E1+ggtitle("(a) Small plasmids"), p.E10, p.E14, nrow=3)
p.big = cowplot::plot_grid(p.FE+ggtitle("(b) Large plasmids"), p.C, nrow=2)
cowplot::plot_grid(p.small, p.big, rel_widths = c(1,1),
                   align='v', axis='tb')
dev.off()

