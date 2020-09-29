#!/n/app/R/3.6.1/bin/Rscript

# Plot correlation between median antisense : sense ratio for coding and tandem antisense

# Written by: K. Lachance
# Date: May 8, 2019

# Use: ./plotCodingTandemCor.R

# Load libraries
library(ggplot2)
library(svglite)

# Read table
coding_dat <- read.table("CodingAntisense/AntisenseSenseRatioMedians.txt", stringsAsFactors = FALSE)
tandem_dat <- read.table("TandemAntisense/AntisenseSenseRatioMedians.txt", stringsAsFactors = FALSE)
labels <- read.table("MutLabels.txt", stringsAsFactors = FALSE)

# Combine tables
dat = cbind(coding_dat, tandem_dat, labels)[,c(1,2,4,6)]
colnames(dat) = c("Mut", "CodingMed", "TandemMed", "Category")

# Calculate correlations
c1 = cor(dat$CodingMed, dat$TandemMed, method = "pearson")
cat(paste("There is a Pearson correlation of ", round(c1, 2), " between median coding sense : antisense ratios and median tandem sense : antisense ratios.\n", sep=""))

# Plot
p1 <- ggplot(dat, aes(x = CodingMed, y = TandemMed, color = Category, label = Mut)) + 
  geom_point() + 
  geom_text() + 
  scale_color_manual(values = c("RNAprocessing" = "#39B54A",
  			        "HistoneVariant" = "#939598",
				"ChromatinModifer" = "#F15A29", 
				"ChromatinRemodeler" = "#1C75BC", 
				"TranscriptionElongation" = "#662D91", 
				"Wildtype" = "black")) + 
  xlab("Median coding ratio") + ylab("Median tandem ratio") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(legend.position = "none") 

ggsave(plot = p1, filename = "CodingTandemCor.svg", width = 5, height = 5)
