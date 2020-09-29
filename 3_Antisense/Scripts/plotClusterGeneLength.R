#!/n/app/R/3.6.1/bin/Rscript

# Plot the distribution od gene lengths for all genes and those in the Rco1 / Set2 cluster

# Written by: K. Lachance
# Date: May 8, 2019

# Use: ./plotClusterGeneLength.R

# Load libraries
library(ggplot2)
library(svglite)

# Read table
cluster_genes <- read.table("Set2Rco1Genes.txt", stringsAsFactors = FALSE)
cluster_genes$Label = "Cluster"
cluster_genes$Length = cluster_genes[,3] - cluster_genes[,2] # End - start

pos_genes <- read.table("../0_Annotations/Genes/genes_polII.noOver.coding_pos.bed", stringsAsFactors = FALSE)
neg_genes <- read.table("../0_Annotations/Genes/genes_polII.noOver.coding_neg.bed", stringsAsFactors = FALSE)
all_genes = rbind(pos_genes, neg_genes)
all_genes$Label = "All"
all_genes$Length = all_genes[,3] - all_genes[,2] # End - start

# Combine tables
dat = rbind(cluster_genes, all_genes)

# Perform T-test to assess statistical difference in distributions
t.test(cluster_genes$Length, all_genes$Length)


# Plot
p1 <- ggplot(dat, aes(x = Label, y = Length/1000, fill = Label)) + 
  geom_boxplot(notch = TRUE, outlier.shape = NA) + 
  scale_fill_manual(values = c("Cluster" = "#507BA7", "All" = "grey50")) + 
  coord_cartesian(ylim = c(0, 5)) + 
  xlab("") + ylab("Gene length (kb)") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(legend.position = "none") 

ggsave(plot = p1, filename = "clusterGeneLengthBoxplot.svg", width = 5, height = 7)
