d = read.csv('../output/targets-genera-n100.csv', header=T)
library(dplyr)
d$length = nchar(d$sequence)
df = d %>% group_by(sequence, genus, length) %>% summarise(n=length(genus))
ggplot(df, aes(genus, n, fill=sequence))+
  geom_bar(stat="identity", position="stack")+
  coord_flip()+
  theme(legend.position ="none")+
  facet_wrap(~length)

# If I could do this with a phylogeny of the genera, 
# that would be really cool. The patterns would poop out
# And it would probably show that we only need to consider
# some R-M systems
df.4 = d[which(d$length==4),] %>% group_by(sequence, genus) %>% summarise(n=length(genus))
ggplot(df.4, aes(genus, n, fill=sequence))+
  geom_bar(stat="identity", position="stack")+
  coord_flip()


df.5 = d[which(d$length==5),] %>% 
  group_by(sequence, genus) %>% 
  summarise(n=length(genus))
ggplot(df.5, aes(genus, n, fill=sequence))+
  geom_bar(stat="identity", position="stack")+
  coord_flip()


df.6 = d[which(d$length==6),] %>% 
  group_by(sequence, genus) %>% 
  summarise(n=length(genus))
ggplot(df.6, aes(genus, n, fill=sequence))+
  geom_bar(stat="identity", position="stack")+
  coord_flip()
