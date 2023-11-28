df = read.csv('../output/PTU-E1-target-counts-unique.csv', header=T)
# Trying to plot core vs. accessory densities

library(ggplot2)
library(tidyr)

df = df[order(df$sequence),]

# # We want to know if sequence has 
# df.seqs = df %>% group_by(sequence, within_range) %>%
#   summarise(n=length(within_range))
# ggplot(df.seqs, aes(sequence, n, group=within_range, fill=within_range))+
#   geom_bar(stat="identity", position="stack")
long_df <- pivot_longer(df, 
                        cols = c(acc_density, core_density),
                        names_to = "variable",
                        values_to = "density")
# deambiguate (normalise by ambiguity)
long_df$density_deambig = long_df$density/long_df$ambiguity
# ggplot(long_df, aes(within_range, density_deambig, colour=variable))+
#   ggbeeswarm::geom_quasirandom(aes(group=variable), position="dodge")+
#   facet_wrap(~length, scales="free_y")


ggplot(long_df, aes(category, density_deambig, colour=variable))+
  geom_boxplot()+
  facet_wrap(~length, scales="free_y")

# PTU-E10
ggplot(long_df[which(long_df$length==6),], aes(category, density_deambig, colour=variable))+
  geom_boxplot(outlier.shape = NA)+
  facet_wrap(~length, scales="free_y")+
  ylim(c(0,0.001))

# PTU-E1
ggplot(long_df[which(long_df$length==6),], aes(category, density_deambig, colour=variable))+
  geom_boxplot(outlier.shape = NA)+
  facet_wrap(~length, scales="free_y")+
  ylim(c(0,0.001))




ggplot(head(long_df[which(long_df$length==6),]), aes(category, density_deambig, colour=variable))+
  geom_point()+
  geom_line(aes(group=sequence))+
  facet_wrap(~sequence, scales="free_y")

ggplot(long_df, aes(category, density, colour=variable))+
  geom_boxplot()+
  facet_wrap(~length, scales="free_y")



# Can calculate difference as well
df$diff = df$core_density-df$acc_density # when positive, greater occurrences in core
ggplot(df, aes(category, diff))+
  geom_boxplot()+
  facet_wrap(~length, scales="free_y")


ggplot(long_df[which(long_df$sequence=="CCWGG"), ], aes(category, density, colour=variable))+
  geom_boxplot()
       

# CCWGG
df = read.csv('../test-CCWGG.csv', header=T)
df$PTU = gsub("-target-.*","", gsub(".*\\/", "", df$file))

long_df <- pivot_longer(df, 
                        cols = c(acc_density, core_density),
                        names_to = "variable",
                        values_to = "density")

ggplot(long_df, aes(variable, density, fill=variable))+
  geom_bar(stat="identity")+
  facet_wrap(~PTU)
