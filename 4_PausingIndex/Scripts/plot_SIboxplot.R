#!/n/app/R/3.6.1/bin/Rscript

# Plot boplot comparing splicing indices between wildtype and cac2 delta

# Written by: K. Lachance
# Date: July 14, 2020

# Use: ./plotBoxplot.R

# Load libraries
#install.packages("ggplot2", repos = "http://cran.us.r-project.org")
library(ggplot2)
library(svglite)

# Read table
cac2_data <- read.table("cac2.junctionCounts.txt", stringsAsFactors = FALSE)
colnames(cac2_data) <- c("Gene", "Intron", "Chr", "Start", "End", "Strand", "Five", "Three", "Spliced")
cac2_data$Mut = "Cac2"

wt_data <- read.table("wt.junctionCounts.txt", stringsAsFactors = FALSE)
colnames(wt_data) <- c("Gene", "Intron", "Chr", "Start", "End", "Strand", "Five", "Three", "Spliced")
wt_data$Mut = "Wt"

# Combine
data = rbind(cac2_data, wt_data)

# Calculate SI
data$SI = (2 * data$Spliced) / (data$Five + data$Three)

# T-test
#t.test(data[data$Mut == "Cac2", 'SI'], data[data$Mut == "Wt", 'SI'], paired = TRUE)

# Wilcoxon rank-sum test
wilcox.test(data[data$Mut == "Cac2", 'SI'], data[data$Mut == "Wt", 'SI'], paired = TRUE)

# Plot 
p <- ggplot(data, aes(x = Mut, y = SI, fill = Mut)) + 
  geom_boxplot() + 
  scale_y_log10() + 
  theme_bw() + 
  scale_fill_manual(values = c("Cac2" = "#1B75BB", "Wt" = "grey50")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") 

ggsave(plot = p, filename = "Cac2SplicingBoxplot.svg", height = 6, width = 3)
