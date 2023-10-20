library(dplyr)
library(ggplot2)
library(ggbeeswarm)

ptu_info = read.csv('../data/top50ptus_info.csv', header=T, stringsAsFactors = F, row.names = 1)
top50_PTUs = read.csv('../top50-ptus-with-core.txt', header=F, stringsAsFactors = F)$V1

gc_df = read.csv('../output/PTU_core_accessory_GC.csv', header=T, stringsAsFactors = F)

gc_df_diff = data.frame(gc_df %>% group_by(PTU) %>%
                          summarise(diff=GC[component=="core"]-GC[component=="accessory"]))
rownames(gc_df_diff) = gc_df_diff$PTU
ptu_info$GC.diff = gc_df_diff[rownames(ptu_info),"diff"]

palindrome_density_df = data.frame(kmer=NA, density=NA, component=NA, PTU=NA)

k = 6
file_suffix = paste0('k', k, '_palindrome_densities.txt')
for (ptu in top50_PTUs){
  components = c("core", "accessory")
  print(ptu)
  for (component in components){
    density_file = paste0('../core-accessory-palindromes/', ptu, '_k', k, '_', component, '_palindrome_densities.txt')
    file.df = read.csv(density_file,sep=' ')
    colnames(file.df) = c("kmer", "density")
    file.df$component = component
    file.df$PTU = ptu
    palindrome_density_df = rbind(palindrome_density_df, file.df)
  }
}
palindrome_density_df = palindrome_density_df[-1,]

# Overall palindrome density?
palindrome_density_overall = palindrome_density_df %>% group_by(PTU, component) %>%
  summarise(density=sum(density))

palindrome_density_diff = palindrome_density_overall %>% group_by(PTU) %>%
  summarise(diff=density[component=="core"]-density[component=="accessory"])

palindrome_density_diff$host.range = ordered(ptu_info[as.character(palindrome_density_diff$PTU), "host.range"], 
                                             levels=c("I", "II", "III", "IV", "V", "VI"))
palindrome_density_diff$host.range.numeric = as.numeric(palindrome_density_diff$host.range)

palindrome_density_diff$size = ptu_info[as.character(palindrome_density_diff$PTU), "median.length"]


p.all = ggplot(palindrome_density_diff, aes(host.range, diff, colour=size>5000))+
  geom_boxplot(outlier.shape = NA, aes(group=host.range), colour="black")+
  ggbeeswarm::geom_quasirandom(size=3)+
  theme_bw()+
  xlab("PTU host range")+
  ylab("palindrome density difference (core-accessory)")+
  theme(axis.text=element_text(size=8, colour="black"), axis.title=element_text(size=12))+
  ggtitle("All plasmids")+
  scale_color_manual(values=c("red", "black"))
ggsave(p.all, file='../output/palindrome-density-top50-plot.pdf')

p.5kb = ggplot(palindrome_density_diff[which(palindrome_density_diff$size>5000),], aes(host.range.numeric, diff))+
  geom_jitter(size=3, width = 0.05, height=0)+
  theme_bw()+
  xlab("PTU host range")+
  stat_smooth(method="lm")+
  ylab("palindrome density difference (core-accessory)")+
  theme(axis.text=element_text(size=8, colour="black"), axis.title=element_text(size=12))+
  ggtitle("Plasmids>5kb")+
  guides(color=FALSE)

# Is it explained by GC content? apparently not
palindrome_density_diff$GC.diff = ptu_info[palindrome_density_diff$PTU, "GC.diff"]
p.GC = ggplot(palindrome_density_diff, aes(GC.diff, diff))+
  geom_point()+
  theme_bw()+
  theme(panel.grid = element_blank())+
  stat_smooth(method="lm")+
  theme(axis.text=element_text(size=8, colour="black"), axis.title=element_text(size=12))+
  xlab("GC difference (core - accessory)")+
  ylab("Palindrome density difference (core - accessory)")+
  ggtitle("No strong GC-associated signal")

p.combined = cowplot::plot_grid(p.all, p.5kb,p.GC, nrow=1,  rel_widths = c(1.2, 1, 1))
ggsave(p.combined, file='../output/palindrome-density-top50-plot.pdf', width=12, height=6)


