df = read.csv('../PTU-Y-test.csv', header=T)
# Trying to plot core vs. accessory densities

library(ggplot2)
library(tidyr)

long_df <- pivot_longer(df, 
                        cols = c(acc_density, core_density),
                        names_to = "variable",
                        values_to = "density")
# deambiguate (normalise by ambiguity)
long_df$density_deambig = long_df$density/long_df$ambiguity
ggplot(long_df, aes(within_range, density_deambig, colour=variable))+
  ggbeeswarm::geom_quasirandom(aes(group=variable), position="dodge")+
  facet_wrap(~length, scales="free_y")
