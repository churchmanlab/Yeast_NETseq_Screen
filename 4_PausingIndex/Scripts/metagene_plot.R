#!/n/app/R/3.6.1/bin/Rscript

# Written by: K. Lachance
# Date: October 30, 2018

# Plot the expanded bed files around loci of interest

# Use: ./metagene_plot.R OverlapMatrix.Trimmed100.pos.txt OverlapMatrix.Trimmed100.min.txt wt TSS

# Import libraries
suppressMessages(library(reshape2))
suppressMessages(library(zoo))

# Get arguments (matrices)
args <- commandArgs(trailingOnly = TRUE)
pos_mat_file <- args[1]
min_mat_file <- args[2]
mut <- args[3]
region <- args[4]

# Load tables
pos_mat <- read.table(pos_mat_file)
min_mat <- read.table(min_mat_file)

# Combine tables
mat <- rbind(pos_mat, min_mat)

# Trim tables
mat2 <- mat[,6:ncol(mat)]

# Normalize each row for total coverage
mat2 = mat2 / rowSums(mat2)
colnames(mat2) <- as.character(1:ncol(mat2))
row.names(mat2) <- mat[,4]

# Melt
melted_mat <- melt(as.matrix(mat2))
colnames(melted_mat) = c("Gene", "Pos", "Value")
melted_mat$Mut = mut
melted_mat$Region = region

write.table(melted_mat, paste("./readyToPlot/", mut, "_combinedMeltedMatrix.bed", sep=""), append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
