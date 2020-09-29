#!/n/app/R/3.6.1/bin/Rscript

# Written by: K. Lachance
# Date: May 8, 2020

# Use: ./plotChIPboxplot.R

# Plot boxplot for ChIP-nexus data for H3/H4 acetylation in gcn5, eaf1, wildtype

# Load libraries
library(ggplot2)
library(svglite)
library(plyr)

# Load data
cat("Loading data...\n")

wt_h3 = read.table("ChIPnexusData/wt_H3K23Ac.Coding.noOver.txt", stringsAsFactors = FALSE)
wt_h3$Label = "wt_h3"
wt_h4 = read.table("ChIPnexusData/wt_H4K12Ac.Coding.noOver.txt", stringsAsFactors = FALSE)
wt_h4$Label = "wt_h4"

eaf1_h3 = read.table("ChIPnexusData/eaf1_H3K23Ac.Coding.noOver.txt", stringsAsFactors = FALSE)
eaf1_h3$Label = "eaf1_h3"
eaf1_h4 = read.table("ChIPnexusData/eaf1_H4K12Ac.Coding.noOver.txt", stringsAsFactors = FALSE)
eaf1_h4$Label = "eaf1_h4"

gcn5_h3 = read.table("ChIPnexusData/gcn5_H3K23Ac.Coding.noOver.txt", stringsAsFactors = FALSE)
gcn5_h3$Label = "gcn5_h3"
gcn5_h4 = read.table("ChIPnexusData/gcn5_H4K12Ac.Coding.noOver.txt", stringsAsFactors = FALSE)
gcn5_h4$Label = "gcn5_h4"

# Combine data frames and add column names
data = rbind(wt_h3, wt_h4, eaf1_h3, eaf1_h4, gcn5_h3, gcn5_h4)
colnames(data) = c("Chr", "Start1", "End1", "Gene", "Skip", "Strand", "Start2", "End2", "Value", "Label")

# Reorder for plotting
data$Label = factor(data$Label, levels = c("wt_h3", "wt_h4", "eaf1_h3", "eaf1_h4", "gcn5_h3", "gcn5_h4"))

# Calculate medians
meds <- ddply(data, .(Label), summarise, med = median(Value))
m = 12.15 # Wildtype H3 median

# Plot
cat("Plotting data...\n")

p <- ggplot(data, aes(x = Label, y = log2(Value), group = Label, fill = Label)) + 
  geom_boxplot(outlier.shape = NA, notch = TRUE) + 
  geom_hline(yintercept = m, linetype = "dashed") + 
  ylim(9, 14) + 
  scale_fill_manual(values = c("wt_h3" = "grey50", "wt_h4" = "grey50", "gcn5_h3" = "#a83e1c", "gcn5_h4" = "#a83e1c", "eaf1_h3" = "#f48b68", "eaf1_h4" = "#f48b68")) + 
  theme_bw() +
  xlab("") + ylab("log2(RPKS histone acetylation)") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(legend.position = "none")

ggsave(filename = "ChIPnexus_acetylation_boxplot.svg", plot = p, width = 10, height = 6)
