#!/n/app/R/3.6.1/bin/Rscript

# Plot correlation of pause strengths (for IDR)

# Written by: K. Lachance
# Date: Feb 5, 2020

# Use: ./plot_pauseScore.R mut

# Load libraries
library(ggplot2)
library(idr)
library(svglite)

# Get input parameters
args <- commandArgs(trailingOnly = TRUE)
MUT <- args[1]

# Get data
df = read.table(paste(MUT, "_OverlapScore.txt", sep=""))
colnames(df) = c("Chr", "Start", "End", "Score1", "Score2", "Strand", "Norm1", "Norm2")

# Calculate correlation
c = round(cor(log10(df$Norm1), log10(df$Norm2), method = "pearson")^2, 2)

# Estimate IDR
dat = cbind(log10(df$Norm1), log10(df$Norm2))
res <- est.IDR(dat, mu=3, sigma=1, rho=.9, p=.5)
dat = cbind(dat, res$idr)

mat = data.frame(dat, stringsAsFactors = FALSE)
colnames(mat) = c("Norm1", "Norm2", "IDR")

# Add T/F for if IDR < 1%
mat$Cat = "Yes"
mat[mat$IDR > 0.01, 'Cat'] = "No"

PASS = df[mat$Cat == "Yes",]

write.table(PASS, paste(MUT, "_PASS_pauseScore.txt", sep=""), quote = F, sep = "\t", row.names = F, col.names = F)

# Plot
# Scatter plot log10 scores
p1 <- ggplot(mat, aes(x = Norm1, y = Norm2, color = Cat)) +
   geom_point(alpha = 0.4) +
   scale_color_manual(values = c("Yes" = "#008A8A", "No" = "#3BBFC3")) +
   theme_bw() +
   xlab("Pause Strength (Rep 1)") + ylab("Pause Strength (Rep 2)") +
   ggtitle(paste(MUT, sep="")) +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
   theme(legend.position = "none")

ggsave(plot = p1, file = paste(MUT, "_IDRscatter.svg", sep=""), height = 5, width = 5)
