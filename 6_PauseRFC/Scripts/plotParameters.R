#!/n/app/R/3.6.1/bin/Rscript

# Plot accuracy of model with various parameters

# Written by: K. Lachance
# Date: May 11, 2020

# Use: ./plotParameters.R mut ncat

# Load libraries
library(ggplot2)
library(svglite)

# Get input parameters
args <- commandArgs(trailingOnly = TRUE)
MUT <- args[1]
nCat <- args[2]

# Load data (generated in O2 in /n/groups/churchman/kcl19/Screen/ANALYSIS/Pausing_THISONE/newThreshold5/IDR/RandomForest/ directory)
data = read.table(paste(, MUT, ".randomForestParameters", nCat, ".txt", header = TRUE, stringsAsFactors = TRUE)

# Report optimal parameters and accuracy
data[which(data$ErrMean == min(data$ErrMean)),]
round(1-min(data$ErrMean), 3)

p <- ggplot(data, aes(x = mtry, y = 1- ErrMean, color = factor(ntree))) + 
  geom_line() + 
  geom_point(shape = 1) + 
  scale_color_manual("Number of Trees", values = c("1000" = "red", "1500" = "orange", "2000" = "green", "2500" = "blue")) + 
  theme_bw() + 
  xlab("Number of variables randomly sampled as candidates at each split") + ylab("Accuracy") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
  
ggsave(filename = paste(MUT, ".optimizedParameters", nCat, ".svg", sep=""), plot = p, height = 3, width = 2)