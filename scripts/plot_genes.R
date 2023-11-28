# Geom gene arrow
library(gggenes)

gene_df = read.csv('~/Desktop/PROKKA_11282023/PROKKA_11282023.gff', 
                   comment.char = "#", sep="\t", header=F)
gene_df = gene_df[gene_df$V3=="CDS",]
gene_df$product = gsub(".*product=", "", gene_df$V9)
gene_df$product[gene_df$product=="hypothetical protein"] = ""
gene_df$length = gene_df$V5-gene_df$V4
pdf('../output/2023-11-28-example-plasmid-genes.pdf', width=20, height=4)
p.gene = ggplot(gene_df, aes(xmin = V4, xmax = V5, y = V7, label=product)) +
  geom_gene_arrow() +
  xlim(c(0,50000))+
  ylab("strand")+
  geom_gene_label(aes(label=product))
dev.off()

# from 2023-11-24-leading-region.R take p.PTU.C
pdf('../output/2023-11-28-example-plasmid.pdf', width=8, height=6)
cowplot::plot_grid(p.PTU.C, p.gene, align='hv', nrow=2, axis='h')
dev.off()
