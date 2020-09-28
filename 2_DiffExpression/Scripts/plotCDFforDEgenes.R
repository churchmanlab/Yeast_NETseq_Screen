#!/n/app/R/3.6.1/bin/Rscript

# Plot CDF for number of strains each DE gene is differentially expressed in

# Written by: K. Lachance
# Date: May 8, 2019

# Use: ./plotCDFforDEgenes.R

# Load libraries
library(ggplot2)
library(gridExtra)

# Read table
dat <- read.table("DEgeneCount.txt", stringsAsFactors = FALSE)
colnames(dat) <- c("nStrains", "Gene")

# Plot
p = ggplot(dat, aes(x = nStrains)) + 
  stat_ecdf(geom = "step") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
ggsave(file="CDFforDEgenes.svg", plot=p, width=8, height=5)

# Calculate and report relevant statistics
nUnique = nrow(dat[dat$nStrains == 1, ])
nTotal = nrow(dat)
cat(paste(round((nUnique / nTotal) * 100, 0), "% of DE genes are differentially expressed in only 1 deletion strain!\n", sep=""))

n90p = quantile(dat$nStrains, 0.9)[[1]]
cat(paste("90% of DE genes are differentially expressed in ", n90p + 1, " strains or fewer!\n", sep=""))