#!/n/app/R/3.6.1/bin/Rscript

# Plot boxplot of pause densities

# Written by: K. Lachance
# Date: June 2, 2020

# Use: ./plot_pauseDensity.R

# Load libraries
library(ggplot2)
library(svglite)

# Read data
data = read.table("ALL.pauseDensity.txt", stringsAsFactors = FALSE)
colnames(data) = c("Chr", "Start", "End", "Gene", "Blank", "Strand", "Npause", "Pdensity", "Mut")

# Calculate median pause density
Pdensity_medians = aggregate(Pdensity ~ Mut, data=data, FUN=median)
colnames(Pdensity_medians) = c("Mut", "Med")

# Order by Comb for plotting
data$Mut <- factor(data$Mut, levels = Pdensity_medians$Mut[order(Pdensity_medians$Med)])

# Color
data$Lab = "a"
data[data$Mut == "wt", 'Lab'] = "b"

# Plot
p <- ggplot(data, aes(x = reorder(Mut, Pdensity, FUN = median), y = Pdensity, group = Mut, color = Lab, fill = Lab)) + 
  geom_boxplot(outlier.shape = NA, coef = 0) + 
  scale_fill_manual(values = c("a" = "white", "b" = "darkcyan")) + 
  scale_color_manual(values = c("a" = "darkcyan", "b" = "black")) + 
  coord_cartesian(ylim = c(0, 150)) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "none")

ggsave(file="pauseDensity.svg", plot=p, width=10, height=6)

