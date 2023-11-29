# Figure 3 plot
# Design: multi-panel. Showing 6-bp palindrome density in big conjugative plasmids 
# with a sliding window.
# With a side panel showing the difference between leading region (5kb) and rest 
# of plasmid.


library(ggplot2)
library(dplyr)

ptu_info = read.csv('../data/top50ptus_info.csv', header=T, row.names = 1)

d = read.csv('../output/231124_leading_region_results.csv', header=T, stringsAsFactors = F)

SUBSET_PTUS = c("PTU-B_O_K_Z",
                    "PTU-C",
                    "PTU-E81",
                    "PTU-E84",
                    "PTU-FE",
                    "PTU-FS",
                    "PTU-HI1A",
                    "PTU-I1",
                    "PTU-I2",
                    "PTU-L_M",
                    "PTU-N1",
                    "PTU-X1",
                    "PTU-X3")

# Nice plot of colours, but shows same as below plot
# p.leading.region.diff = ggplot(d[which(d$PTU %in% SUBSET_PTUS),], aes(leading_region_density, rest_plasmid_density))+
#   geom_point(aes(colour=PTU))+
#   xlab("leading region")+
#   ylab("rest")+
#   coord_cartesian(xlim = c(0,0.02), ylim=c(0,0.02))+
#   geom_abline(slope=1)+
#   theme_bw()+
#   theme(panel.grid = element_blank(),
#         axis.text=element_text(colour="black"))+
#   xlab("Leading region (5kb)")+
#   ylab("Rest of plasmid")+
#   ggtitle("6-bp palindrome density")

d$host.range = ptu_info[as.character(d$PTU),"host.range"]

d$diff = d$leading_region_density-d$rest_plasmid_density
d$PTU = ordered(d$PTU, 
                levels=rev(unique(d$PTU[order(d$diff)])))
p.leading.region.diff = ggplot(d[which(d$PTU %in% SUBSET_PTUS),], aes(PTU, diff))+
  geom_hline(yintercept = 0, linetype='dashed')+
  theme_bw()+
  ggbeeswarm::geom_quasirandom()+
  stat_summary(colour="red", fun="median")+
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  facet_wrap(~host.range, scales="free_x", nrow=2)+
  theme(panel.grid = element_blank(),
        axis.text=element_text(colour="black"))+
  xlab("")+
  ylab("6-bp palindrome density in leading region (5kb)\nvs. density in rest of plasmid")
  


# Leading region sliding window
d = read.csv('../output/231124_leading_region_results_window_v2.csv', header=F, stringsAsFactors = F)
d.reduced = d[which((d$V3 %% 500)==0), ]

colnames(d.reduced) = c("PTU", "plasmid", "x", "density", "GC_n", "expected_density")

d.reduced$diff_density = d.reduced$density-d.reduced$expected_density

medians = d.reduced %>% group_by(PTU, plasmid) %>%
  summarise(median=median(density)) %>% group_by(PTU) %>%
  summarise(median=median(median))
# p = ggplot(d.reduced,
#            aes(x, density*5000))+
#   geom_line(aes(group=plasmid), alpha=0.2)+
#   facet_wrap(~PTU, scales="free")

# for (i in 1:nrow(medians)) {
#   p <- p + geom_hline(data=medians[i,], aes(yintercept=median*5000), color="red")
# }

# One line per plasmid
d.reduced.median = d.reduced %>% group_by(PTU, x) %>% 
  summarise(median=median(diff_density),
            n=length(diff_density), 
            median.original=median(density),
            GC=median(GC_n))


d.reduced.median.plot = d.reduced.median[which(d.reduced.median$PTU %in% SUBSET_PTUS),]

# Maximum counts
d.reduced.median.plot = d.reduced.median.plot %>% 
  group_by(PTU) %>%
  mutate(max_n = max(n, na.rm = TRUE),  # Calculate max 'n' for each 'PTU'
         normalized_n = n / max_n * 100) %>%  # Normalize 'n' by this max value
  ungroup() %>%  # Remove the grouping
  select(-max_n)  # Optionally remove the max_n column if not needed


p.lines = ggplot(d.reduced.median.plot, aes(x, median.original))+
  geom_hline(yintercept = 0, linetype='dashed')+
  theme_bw()+
  geom_line(aes(group=PTU, alpha=normalized_n))+
  facet_wrap(~PTU, scales="free_x", nrow=3)+
  ylab("6-mer palindrome density (5kb window)")+
  labs(alpha="% of plasmids within PTU")+
  theme(panel.grid=element_blank())+  
  scale_x_continuous(labels = function(x) x / 1000) +
  xlab("Position on plasmid (kb)")+
  theme(legend.position = c(0.8, 0.2))
  



pdf(paste0('../figures/', Sys.Date(), '-Figure-3.pdf'), width=15, height=8)
cowplot::plot_grid(p.lines+ggtitle("(a)"), 
                   p.leading.region.diff+ggtitle("(b)"), rel_widths = c(1,0.5),
                   nrow=1, align='v', axis='bt')
dev.off()

