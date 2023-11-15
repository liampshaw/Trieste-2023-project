d = read.csv('../output/targets-genera-n100.csv', header=T)
library(dplyr)
library(ggplot2)
d$length = nchar(d$sequence)
df = d %>% group_by(sequence, genus, length) %>% summarise(n=length(genus))
ggplot(df, aes(genus, n, fill=sequence))+
  geom_bar(stat="identity", position="stack")+
  coord_flip()+
  theme(legend.position ="none")+
  facet_wrap(~length)


library(ape)
library(phangorn)
library(ggtree)
t = read.tree('../data/PTU-genera-16S.tre')
t = midpoint(t)
plot(t)

library(AMR)
df$class = mo_class(df$genus)

tax.df = data.frame(genus=unique(df$genus),
                    family=mo_family(unique(df$genus)),
                    order=mo_order(unique(df$genus)),
                    class=mo_class(unique(df$genus)),
                    phylum=mo_phylum(unique(df$genus)))

tax.df = rbind(tax.df,
               c("Chlamydia",
                 mo_family("Chlamydia"),
                 mo_order("Chlamydia"),
                 mo_class("Chlamydia"),
                 mo_phylum("Chlamydia")))
t.outed = root(t, "Chlamydia")
p = ggtree(t.outed)
p.tree = p %<+% tax.df+
  geom_tiplab(size=1, align = TRUE)

p.tree.plot.class = p %<+% tax.df+
  geom_tiplab(size=3, align = TRUE)+
  geom_tippoint(aes(colour=class))
p.tree.plot.order = p %<+% tax.df+
  geom_tiplab(size=3, align = TRUE)+
  geom_tippoint(aes(colour=order))
pdf('../output/2023-11-14-genera-phylogeny-class.pdf')
p.tree.plot.class
dev.off()

pdf('../output/2023-11-14-genera-phylogeny-order.pdf')
p.tree.plot.order
dev.off()

GENUS.LEVELS = rev(c("Escherichia", "Shigella", "Salmonella", "Cronobacter", "Pantoea", "Erwinia", "Pluralibacter", "Kluyvera", "Klebsiella", "Citrobacter", "Enterobacter", "Tatumella", "Serratia", "Raoultella", "Pectobacterium", "Providencia", "Proteus", "Morganella", "Edwardsiella", "Yersinia", "Aeromonas", "Vibrio", "Photobacterium", "Shewanella", "Acinetobacter", "Burkholderia", "Staphylococcus", "Macrococcus", "Bacillus", "Streptococcus", "Enterococcus", "Corynebacterium", "Rhodopseudomonas", "Campylobacter", "Chlamydia"))
df$genus = ordered(df$genus,
                   levels=GENUS.LEVELS)
ggplot(na.omit(df), aes(genus, n, fill=sequence))+
  geom_bar(stat="identity", position="stack")+
  coord_flip()+
  theme(legend.position ="none")+
  facet_wrap(~length)

# If I could do this with a phylogeny of the genera, 
# that would be really cool. The patterns would poop out
# And it would probably show that we only need to consider
# some R-M systems
df.4 = df[which(df$length==4),]
df.4$genus = ordered(df.4$genus,
                     levels=GENUS.LEVELS)
p.4 = ggplot(df.4, aes(genus, n, fill=sequence))+
  geom_bar(stat="identity", position="stack")+
  coord_flip()+
  scale_x_discrete(drop=FALSE)+
  xlab("")+
  theme_bw()

df.5 = df[which(df$length==5),]
df.5$genus = ordered(df.5$genus,
                     levels=GENUS.LEVELS)
p.5 = ggplot(na.omit(df.5), aes(genus, n, fill=sequence))+
  geom_bar(stat="identity", position="stack")+
  coord_flip()+
  scale_x_discrete(drop=FALSE)+
  xlab("")+
  theme_bw()


df.6 = df[which(df$length==6),]
df.6$genus = ordered(df.6$genus,
                     levels=GENUS.LEVELS)
p.6 = ggplot(na.omit(df.6), aes(genus, n, fill=sequence))+
  geom_bar(stat="identity", position="stack")+
  coord_flip()+
  scale_x_discrete(drop=FALSE)+
  xlab("")+
  theme_bw()


pdf('../output/2023-11-14-k4-target-db.pdf')
cowplot::plot_grid(p.tree, p.4, 
                   nrow=1, align='h',
                   rel_widths = c(1, 2))
dev.off()

pdf('../output/2023-11-14-k5-target-db.pdf')
cowplot::plot_grid(p.tree, p.5, 
                   nrow=1, align='h',
                   rel_widths = c(1, 4))
dev.off()


pdf('../output/2023-11-14-k6-target-db.pdf')
cowplot::plot_grid(p.tree, p.6, 
                   nrow=1, align='h',
                   rel_widths = c(1, 4))
dev.off()




# Investigating escherichia vs. shigella
# 14/11/2023
escherichia.ptus = unique(d$PTU[which(d$Genus=="Escherichia")])
shigella.ptus = unique(d$PTU[which(d$Genus=="Shigella")])

overlap.ptus = escherichia.ptus[which(escherichia.ptus %in% shigella.ptus)]
d.overlap = d[which(d$PTU %in% overlap.ptus & d$Genus %in% c("Escherichia", "Shigella")),]
d.overlap.summary = d.overlap %>% group_by(PTU) %>%
  summarise(size=median(Size), n.escherichia=length(Genus[Genus=="Escherichia"]),
            n.shigella=length(Genus[Genus=="Shigella"]))
