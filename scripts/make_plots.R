
PTU.LEVELS = c("PTU-FS", "PTU-I1" , "PTU-C", "PTU-E1", "PTU-Q1")
ptus = c("PTU-FS", "PTU-I1", "PTU-C", "PTU-E1", "PTU-Q1")

all_plasmids = data.frame(kmer=NA, density=NA, plasmid=NA, PTU=NA)

k = 6
pattern_to_search = paste0('k', k, '_palindrome_densities.txt')
for (ptu in ptus){
  density_files = list.files(pattern=pattern_to_search,ptu)
  for (file in density_files){
    file.df = read.csv(paste0(ptu, "/",file),sep=' ')
    colnames(file.df) = c("kmer", "density")
    plasmid_name = gsub("_k.*", "", file)
    file.df$plasmid = plasmid_name
    file.df$PTU = ptu
    all_plasmids = rbind(all_plasmids, file.df)
  }
}

all_plasmids = all_plasmids[-1,]


all_plasmids_k4 = data.frame(kmer=NA, density=NA, plasmid=NA, PTU=NA)
k = 4
pattern_to_search = paste0('k', k, '_palindrome_densities.txt')
for (ptu in ptus){
  density_files = list.files(pattern=pattern_to_search,ptu)
  for (file in density_files){
    file.df = read.csv(paste0(ptu, "/",file),sep=' ')
    colnames(file.df) = c("kmer", "density")
    plasmid_name = gsub("_k.*", "", file)
    file.df$plasmid = plasmid_name
    file.df$PTU = ptu
    all_plasmids_k4 = rbind(all_plasmids_k4, file.df)
  }
}

all_plasmids_k4 = all_plasmids_k4[-1,]

library(ggplot2)
library(ggbeeswarm)
library(dplyr)



ggplot(all_plasmids, aes(PTU, density))+
  ggbeeswarm::geom_quasirandom()+
  facet_wrap(~kmer)+
  theme(axis.text=element_text(angle=45, hjust=1))+
  theme_bw()


# Overall palindrome density?
all_plasmids_density_overall = all_plasmids %>% group_by(plasmid, PTU) %>%
  summarise(density=sum(density))


all_plasmids_density_overall$PTU = ordered(all_plasmids_density_overall$PTU,
                                           levels=PTU.LEVELS)
p.k6 = ggplot(all_plasmids_density_overall, aes(PTU, density))+
  ggbeeswarm::geom_quasirandom()+
  theme(axis.text.x=element_text(angle=45, hjust=1, size=18),
        axis.text.y=element_text(size=18))+
  theme_bw()+
  theme(panel.grid=element_blank())+
  ylab("Overall palindrome density")+
  xlab("")+
  ggtitle("k=6")


# Overall palindrome density k=4
all_plasmids_density_overall_k4 = all_plasmids_k4 %>% group_by(plasmid, PTU) %>%
  summarise(density=sum(density))

all_plasmids_density_overall_k4$PTU = ordered(all_plasmids_density_overall_k4$PTU,
                                           levels=PTU.LEVELS)

p.k4 = ggplot(all_plasmids_density_overall_k4, aes(PTU, density))+
  ggbeeswarm::geom_quasirandom()+
  theme(axis.text.x=element_text(angle=45, hjust=1, size=18),
        axis.text.y=element_text(size=18))+
  theme_bw()+
  ylab("Overall palindrome density")+
  xlab("")+
  theme(panel.grid=element_blank())+
  ggtitle("k=4")

ggsave(p.k4, width=6, height=4, file='plots/overall-densities-k4.pdf')
ggsave(p.k6, width=6, height=4, file='plots/overall-densities-k6.pdf')


# CORE vs. ACCESSORY
k4_gene_densities = data.frame(kmer=NA, density=NA, component=NA, PTU=NA, k=NA)
k6_gene_densities =  data.frame(kmer=NA, density=NA,  component=NA, PTU=NA, k=NA)
for (ptu in ptus){
  for (component in c("core", "accessory")){
    file.df = read.csv(paste0(ptu, "/",'k',4, '_', component, '_palindrome_densities.txt' ),sep=' ')
    colnames(file.df) = c("kmer", "density")
    plasmid_name = gsub("_k.*", "", file)
    file.df$k = 4
    file.df$PTU = ptu
    file.df$component = component
    k4_gene_densities = rbind(k4_gene_densities, file.df)
    file.df = read.csv(paste0(ptu, "/",'k',6, '_', component, '_palindrome_densities.txt' ),sep=' ')
    colnames(file.df) = c("kmer", "density")
    plasmid_name = gsub("_k.*", "", file)
    file.df$k = 6
    file.df$PTU = ptu
    file.df$component = component
    
    k6_gene_densities = rbind(k6_gene_densities, file.df)
  }
    
}
# drop NA rows
k4_gene_densities   = k4_gene_densities[-1,]
k6_gene_densities   = k6_gene_densities[-1,]

k4_gene_densities$PTU = ordered(k4_gene_densities$PTU,
                                levels=PTU.LEVELS)
k6_gene_densities$PTU = ordered(k6_gene_densities$PTU,
                                levels=PTU.LEVELS)
k4_gene_densities_overall = k4_gene_densities %>% group_by(PTU, component) %>%
  summarise(density=sum(density))
p.k4.genes = ggplot(k4_gene_densities_overall, aes(group=component, PTU, density, fill=component))+
  geom_bar(stat="identity", position="dodge")+
  theme_bw()+
  theme(panel.grid=element_blank())+
  scale_fill_manual(values=c("#1b9e77", "#d95f02"))+
  ylab("Overall palindrome density")+
  xlab("")+
  theme(axis.text=element_text( size=12),
        axis.text.x=element_text(size=12, angle=45, hjust=1),
        axis.title=element_text(size=18),
        strip.text=element_text(size=18),
        title=element_text(size=18))+
  ggtitle("k=4")+
  facet_wrap(~component)+
  theme(legend.position = "none")


k6_gene_densities_overall = k6_gene_densities %>% group_by(PTU, component) %>%
  summarise(density=sum(density))
p.k6.genes = ggplot(k6_gene_densities_overall, aes(group=component, PTU, density, fill=component))+
  geom_bar(stat="identity", position="dodge")+
  theme_bw()+
  theme(panel.grid=element_blank())+
  scale_fill_manual(values=c("#1b9e77", "#d95f02"))+
  ylab("Overall palindrome density")+
  xlab("")+
  theme(axis.text=element_text( size=12),
        axis.text.x=element_text(size=12, angle=45, hjust=1),
        axis.title=element_text(size=18),
        strip.text=element_text(size=18),
        title=element_text(size=18))+
  ggtitle("k=6")+
  facet_wrap(~component)+
  theme(legend.position = "none")
ggsave(p.k4.genes, width=6, height=4, file='plots/gene-densities-k4.pdf')
ggsave(p.k6.genes, width=6, height=4, file='plots/gene-densities-k6.pdf')


k6_gene_densities_diff = k6_gene_densities %>% group_by(PTU, component) %>%
  summarise(density=sum(density)) %>% 
  group_by(PTU) %>%
  summarise(diff=density[component=="core"]-density[component=="accessory"])

k4_gene_densities_diff = k4_gene_densities %>% group_by(PTU, component) %>%
  summarise(density=sum(density)) %>% 
  group_by(PTU) %>%
  summarise(diff=density[component=="core"]-density[component=="accessory"])
