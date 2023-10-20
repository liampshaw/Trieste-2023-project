library(dplyr)
library(ggplot2)
library(ggbeeswarm)

ptu_info = read.csv('../data/top50ptus_info.csv', header=T, stringsAsFactors = F, row.names = 1)
top50_PTUs = read.csv('../top50-ptus-with-core.txt', header=F, stringsAsFactors = F)$V1
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


ggplot(palindrome_density_diff, aes(host.range, diff))+
  geom_boxplot(outlier.shape = NA)+
  ggbeeswarm::geom_quasirandom()+
  theme_bw()+
  xlab("PTU host range")+
  ylab("palindrome density difference (core-accessory)")+
  theme(axis.text=element_text(size=12, colour="black"), axis.title=element_text(size=18))+
  ggtitle("All plasmids")
ggplot(palindrome_density_diff[which(palindrome_density_diff$size>5000),], aes(host.range, diff))+
  geom_boxplot(outlier.shape = NA)+
  ggbeeswarm::geom_quasirandom()+
  theme_bw()+
  xlab("PTU host range")+
  ylab("palindrome density difference (core-accessory)")+
  theme(axis.text=element_text(size=12, colour="black"), axis.title=element_text(size=18))+
  ggtitle("Median size > 5kb")
ggplot(palindrome_density_diff[which(palindrome_density_diff$size>20000),], aes(host.range, diff))+
  geom_boxplot(outlier.shape = NA)+
  ggbeeswarm::geom_quasirandom()+
  theme_bw()+
  xlab("PTU host range")+
  ylab("palindrome density difference (core-accessory)")+
  theme(axis.text=element_text(size=12, colour="black"), axis.title=element_text(size=18))+
  ggtitle("Median size > 20kb")

summary(lm(diff ~ host.range.numeric:log10(size), data=palindrome_density_diff))

ggplot(palindrome_density_overall[which(palindrome_density_overall$PTU=="PTU-C"),], aes(component, density))+
  geom_bar(stat="identity")
