library(ggplot2)

ptu_info = read.csv('../data/top50ptus_info.csv', header=T, row.names = 1)

d = read.csv('../output/231124_leading_region_results.csv', header=T, stringsAsFactors = F)


ggplot(d, aes(leading_region_density, rest_plasmid_density))+
  geom_point(aes(colour=PTU))+
  xlab("leading region")+
  ylab("rest")+
  coord_cartesian(xlim = c(0,0.02), ylim=c(0,0.02))+
  geom_abline(slope=1)

d$host.range = ptu_info[as.character(d$PTU),"host.range"]

d$diff = d$leading_region_density-d$rest_plasmid_density
d$PTU = ordered(d$PTU, 
                levels=rev(unique(d$PTU[order(d$diff)])))
ggplot(d, aes(PTU, diff))+
  geom_hline(yintercept = 0, linetype='dashed')+
  theme_bw()+
  geom_boxplot(aes(colour=PTU))+
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  facet_wrap(~host.range, scales="free_x", nrow=1)

ggplot(d, aes(plasmid_length, diff))+
  geom_point()+
  facet_wrap(~host.range)


ggplot(d, aes(leading_region_GC, rest_of_plasmid_GC))+
  geom_point()


ggplot(d, aes(leading_region_GC-rest_of_plasmid_GC,diff))+
  geom_point()+
  ylab("diff between palindrome density")+
  facet_wrap(~plasmid)


# Leading region sliding window
d = read.csv('../output/231124_leading_region_results_window_v2.csv', header=F, stringsAsFactors = F)
d.reduced = d[which((d$V3 %% 500)==0), ]

colnames(d.reduced) = c("PTU", "plasmid", "x", "density", "GC_n", "expected_density")
library(ggplot2)
library(dplyr)

medians = d.reduced %>% group_by(PTU, plasmid) %>%
  summarise(median=median(density)) %>% group_by(PTU) %>%
  summarise(median=median(median))
p = ggplot(d.reduced,
       aes(x, density*5000))+
  geom_line(aes(group=plasmid), alpha=0.2)+
  facet_wrap(~PTU, scales="free")

for (i in 1:nrow(medians)) {
  p <- p + geom_hline(data=medians[i,], aes(yintercept=median*5000), color="red")
}

pdf("../output/2023-11-24-sliding-window.pdf", height=20, width=20)
p
dev.off()

d.reduced$diff_density = d.reduced$density-d.reduced$expected_density
ggplot(d.reduced[which(d.reduced$PTU=="PTU-I1"),], aes(x, diff_density))+
  geom_hline(yintercept = 0, linetype='dashed')+
  theme_bw()+
  geom_line(aes(group=plasmid), alpha=0.1)

# One line per plasmid?
d.reduced.median = d.reduced %>% group_by(PTU, x) %>% 
  summarise(median=median(diff_density))
ggplot(d.reduced.median, aes(x, median))+
  geom_hline(yintercept = 0, linetype='dashed')+
  theme_bw()+
  geom_line(aes(group=PTU), alpha=1)+
  facet_wrap(~PTU, scales="free_x")
