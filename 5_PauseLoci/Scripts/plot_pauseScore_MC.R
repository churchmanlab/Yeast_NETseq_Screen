#!/n/app/R/4.0.1/bin/Rscript

# Plot correlation of pause strengths (for IDR)

# Written by: K. Lachance
# Date: Feb 5, 2020
# M.Couvillion: 6/7/21: update IDR_reproducibility.txt to contain number passing IDR

# Use: ./plot_pauseScore.R mut covT

# Load libraries
library(ggplot2)
library(idr)
library(svglite)
library(data.table)

# Get input parameters
args <- commandArgs(trailingOnly = TRUE)
MUT <- args[1]
covT <- args[2]

# Get data
df = read.table(paste0(MUT, "_OverlapScore_cov",covT,".txt"))
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

write.table(PASS, paste(MUT, "_PASS_pauseScore_cov",covT,".txt", sep=""), quote = F, sep = "\t", row.names = F, col.names = F)

# Add number passing IDR to IDR_reproducibility.txt
IDRpass = nrow(PASS)

IDRrep_DT = data.table(read.table(paste0("IDR_reproducibility_cov", covT, "_tmp.txt")))
line = IDRrep_DT[IDRrep_DT$V1 == MUT]
line[, IDR := IDRpass]
write.table(line, paste0("IDR_reproducibility_cov", covT, ".txt"), append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# Plot
# Scatter plot log10 scores
p1 <- ggplot(mat, aes(x = Norm1, y = Norm2, color = Cat)) +
   geom_point(alpha = 0.4) +
   scale_color_manual(values = c("Yes" = "#008A8A", "No" = "#3BBFC3")) +
   coord_cartesian(ylim = c(-1, 2.2)) + 
   theme_bw() +
   xlab("Pause Strength (Rep 1)") + ylab("Pause Strength (Rep 2)") +
   ggtitle(paste(MUT, sep="")) +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
   theme(legend.position = "none")

ggsave(plot = p1, file = paste0(MUT, "_IDRscatter_cov", covT,".svg"), height = 3, width = 3)
