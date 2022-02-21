#!/n/app/R/4.0.1/bin/Rscript

# Plot trends in pauses across mutants

# Written by: K. Lachance
# Date: May 8, 2019

# Use: ../Scripts/plotNumDEgenes_subgenic.R

# Load libraries
library(ggplot2)
library(reshape2)
library(gridExtra)

# Read table
dat <- read.table("Screen_sub_DE.txt", header = FALSE, stringsAsFactors = FALSE)
dat2 <- read.table("Screen_DE.txt", header = FALSE, stringsAsFactors = FALSE)
colnames(dat) <- c("Mut", "Up", "Down", "Total")
colnames(dat2) <- c("Mut", "Up", "Down", "Total")

# Reorder for plotting
dat$Mut <- reorder(dat2$Mut, dat2$Total)

# Melt data
m_dat = melt(dat, id.vars = "Mut", measure.vars = c("Up", "Down"))

# Plot
p1 <- ggplot(m_dat, aes(x = Mut, y = value, fill = variable)) + 
    geom_bar(stat = "identity") + 
    scale_fill_manual(values=c("Up" = "dodgerblue3", "Down" = "coral3")) + 
    ylab("Number of DE Genes") + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    theme(legend.position="none")

ggsave(file="NumDEgenes_subgenic.pdf", plot=p1, width=8, height=5)

