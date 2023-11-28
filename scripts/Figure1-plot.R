# Figure 1 plot
# Design: three panels. From left to right:
# - phylogeny
# - plasmid prevalence heatmap per PTU (number seen per species)
# - R-M system distribution (6-mers only)

library(ape)
library(phangorn)
library(ggplot2)
library(ggtree)
library(dplyr)
library(tidyr)
library(scales) 
library(aplot)

# Custom colour function for heatmap
black_color_scale <- scale_fill_gradientn(
  colours = c("white", "darkgrey", "black"),
  values = scales::rescale(c(0, 0.01, 1), to = c(0, 1)),
  na.value = "white"
)
red_color_scale <- scale_fill_gradientn(
  colours = c("white", "pink", "red"),
  values = scales::rescale(c(0, 0.01, 1), to = c(0, 1)),
  na.value = "white"
)

# Read in tree
tree = read.tree('../data/PTU-genera-16S.tre')
tree = midpoint(tree)

  
# Read in plasmids
plasmid_df = read.csv('../data/table_50PTUs_RS200.csv',
                      header=T,
                      sep='\t',
                      row.names = 1)
# Omit 4 no genus plasmids
# length(which(plasmid_df$Genus=="-")) 
# Omit 9 Leclercia plasmids
# length(which(plasmid_df$Genus=="Leclercia")) 
# Omit 2 Phytobacter plasmids
# length(which(plasmid_df$Genus=="Phytobacter")) 
plasmid_df = plasmid_df[which(plasmid_df$Genus!="-" &
                                plasmid_df$Genus!="Leclercia" & 
                                plasmid_df$Genus!="Phytobacter"),]
plasmids_by_genera = plasmid_df %>% group_by(Genus, PTU) %>%
  summarise(n=length(Genus))

# Sort PTUs by number of distinct genera
ptu_counts <- plasmids_by_genera %>%
  group_by(PTU) %>%
  summarise(Genus_Count = n_distinct(Genus)) %>%
  arrange(desc(Genus_Count))

# Combine this as a heatmap with tree - needs to be data.frame
ordered_data <- data.frame(plasmids_by_genera %>%
  pivot_wider(names_from = PTU, values_from = n, values_fill = list(n = 0)) %>%
  select(c("Genus", one_of(ptu_counts$PTU))))

# Make sure row names are the genus names for gheatmap
rownames(ordered_data) <- ordered_data$Genus
ordered_data$Genus <- NULL
colnames(ordered_data) = ptu_counts$PTU

# Just presence/absence

  

color_vector <- unlist(lapply(heatmap_data, 
                              function(column) sapply(column, 
                                                      get_color)))
get_color(2)

# Plot
# Rename PTUs to remove 'PTU-' for ease of plotting
tree.ptus = tree 
tree.ptus$tip.label = gsub("PTU-", "", tree$tip.label)
# Plot phylogeny 
p.tree = ggtree(tree.ptus)+
  geom_tiplab(size=3)

# Ordered data with new column names
colnames(ordered_data) = gsub("PTU-", "", colnames(ordered_data))

# And heatmap of PTUs
p.2 = gheatmap(p.tree, ordered_data, 
          offset = 0.1, 
          width = 2.5, 
          colnames_angle = 45,
          hjust=1,
          font.size = 3
          )+
  black_color_scale+
  labs(fill="Plasmids")
p.2 = p.2 + ggtree::vexpand(.05, -1)

ggsave(p.2, filename='~/Desktop/test.pdf', width=15, height = 7)
  
# Read in RM target data
rm_target_df = read.csv('../output/targets-genera-n100.csv', header=T)
# Omit Leclercia and no genus
rm_target_df = rm_target_df[which(rm_target_df$genus!="-" & 
                                    rm_target_df$genus!="Leclercia" & 
                                    rm_target_df$genus!="Phytobacter"),]

rm_target_df$length = nchar(rm_target_df$sequence)
rm_target_df_6 = rm_target_df %>% 
  filter(length==6) %>%
  group_by(sequence, genus, length) %>% 
  summarise(n=length(genus))
rm_target_df_6$length = NULL
colnames(rm_target_df_6) = c("Sequence", "Genus", "n")
# Sort by number of genera RM systems are seen in
rm_counts = rm_target_df_6 %>%
  group_by(Sequence) %>%
  summarise(Genus_Count = n_distinct(Genus)) %>%
  arrange(desc(Genus_Count))

rm_heatmap_data = data.frame(rm_target_df_6 %>%
             pivot_wider(names_from = Sequence, values_from = n, values_fill = list(n = 0)) %>%
             select(c("Genus", one_of(rm_counts$Sequence))))
rownames(rm_heatmap_data) = rm_heatmap_data$Genus
rm_heatmap_data$Genus = NULL

rm_target_df_6$Genus = ordered(rm_target_df_6$Genus,
                               levels=)
plasmids_by_genera$PTU = ordered(plasmids_by_genera$PTU,
                                 levels=ptu_counts$PTU)
p.ptu = ggplot(plasmids_by_genera,aes(PTU,Genus, fill=n) )+
  geom_tile()+
  theme_bw()+
  theme(panel.grid = element_blank())+
  black_color_scale+
  theme(axis.text=element_text(colour="black", size=6,
                               angle=45, hjust=1))+
  labs(fill="No. plasmids")+
  ylab("")+
  theme(axis.text.y=element_blank())

rm_target_df_6$Sequence = ordered(rm_target_df_6$Sequence,
                                 levels=rm_counts$Sequence)
p.rm = ggplot(rm_target_df_6, aes(Sequence, Genus, fill=n))+
  geom_tile()+
  red_color_scale+
  theme_bw()+
  theme(panel.grid=element_blank())+
  theme(axis.text=element_text(colour="black", 
                               angle=45, 
                               hjust=1,
                              size=6),
        axis.text.y=element_blank())+
  labs(fill="% prevalence of systems")+
  ylab("")+
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "points"),  # Remove plot margin
    panel.spacing = unit(0, "lines"),  # Remove spacing between panels
  )+
  xlab("R-M target")

p.combined = p.ptu %>% 
  insert_left(p.tree+ggtree::hexpand(.7, 1), width=0.5) %>% 
  insert_right(p.rm, width=1.1)
pdf('../figures/2023-11-28-Figure-1.pdf', width=15, height=6)
p.combined
dev.off()

