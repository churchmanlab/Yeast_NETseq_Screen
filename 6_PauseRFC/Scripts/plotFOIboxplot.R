#!/n/app/R/3.6.1/bin/Rscript

# Plot distributions of features used to discriminate between real and shuffled pauses

# Written by: K. Lachance
# Date: May 11, 2020

# Use: ./plotFeatures.R mut

# Load libraries
cat("Loading libraries...\n")
library(ggplot2)
library(svglite)

# Get input parameters
args <- commandArgs(trailingOnly = TRUE)
MUT <- args[1]

### JUST ONE FEATURE OF INTEREST
FOI = args[2]

# Load data
cat("Formatting data...\n")
data_all = read.table(paste(MUT, ".IDRrep.PAUSESwithFEATURES.bed", sep=""), header = TRUE, stringsAsFactors = TRUE)

t.test(x = data_all[data_all$Real == 'Real', FOI], y = data_all[data_all$Real == 'Shuffle', FOI], alternative = "two.sided")

data_all = data_all[data_all[,FOI] > 0,]

p1 <- ggplot(data_all, aes(x = Real, y = get(FOI)/1000, group = Real, fill = Real)) +
   geom_boxplot(notch = TRUE, outlier.shape = NA, coef = 0) + 
   scale_fill_manual(values = c("navy", "grey50")) +
   theme_bw() + 
   theme(legend.position = "none") +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
   ggtitle(paste(FOI)) + theme(plot.title = element_text(hjust = 0.5))

ggsave(filename = paste(MUT, "_", FOI, ".featureDistribution.svg", sep=""), plot = p1, height = 3, width = 2)
