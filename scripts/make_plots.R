
ptus = c("PTU-I1", "PTU-C", "PTU-FS", "PTU-E1", "PTU-Q1")

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
                                           levels=c("PTU-I1", "PTU-C", "PTU-FS", "PTU-E1", "PTU-Q1"))
p.k6 = ggplot(all_plasmids_density_overall, aes(PTU, density))+
  ggbeeswarm::geom_quasirandom()+
  theme(axis.text=element_text(angle=45, hjust=1))+
  theme_bw()+
  ylab("Overall palindrome density")+
  xlab("")+
  ggtitle("k=6")


# Overall palindrome density k=4
all_plasmids_density_overall_k4 = all_plasmids_k4 %>% group_by(plasmid, PTU) %>%
  summarise(density=sum(density))

all_plasmids_density_overall_k4$PTU = ordered(all_plasmids_density_overall_k4$PTU,
                                           levels=c("PTU-I1", "PTU-C", "PTU-FS", "PTU-E1", "PTU-Q1"))

p.k4 = ggplot(all_plasmids_density_overall_k4, aes(PTU, density))+
  ggbeeswarm::geom_quasirandom()+
  theme(axis.text=element_text(angle=45, hjust=1))+
  theme_bw()+
  ylab("Overall palindrome density")+
  xlab("")+
  ggtitle("k=4")

ggsave(p.k4, file='plots/overall-densities-k4.pdf')
ggsave(p.k6, file='plots/overall-densities-k6.pdf')

  
