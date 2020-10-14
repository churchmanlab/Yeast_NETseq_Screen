#!/n/app/R/3.6.1/bin/Rscript

# Plot correlation between AUC and number of pauses

# Written by: K. Lachance
# Date: May 8, 2019

# Use: ./plotAUCnPauseCor.R

# Load libraries
library(ggplot2)
library(svglite)

# Read table
dat <- read.table("AUC_nPause.txt", stringsAsFactors = FALSE)
labels <- read.table("MutLabels.txt", stringsAsFactors = FALSE)

# Combine tables
dat = cbind(dat, labels)[,c(1,2,4)]
colnames(dat) = c("Mut", "AUC", "nPause", "Category")

# Calculate correlations
c1 = cor(dat$AUC, dat$nPause, method = "pearson")
cat(paste("There is a Pearson correlation of ", round(c1, 2), " between AUC values and number of pauses in each deletion strain.\n", sep=""))
cor.test(dat$AUC, dat$nPause)

# Plot
p1 <- ggplot(dat, aes(x = AUC, y = nPause / 1000, color = Category, label = Mut)) + 
  geom_point() + 
  geom_text() + 
  scale_color_manual(values = c("RNAprocessing" = "#39B54A",
  			        "HistoneVariant" = "#939598",
				"ChromatinModifer" = "#F15A29", 
				"ChromatinRemodeler" = "#1C75BC", 
				"TranscriptionElongation" = "#662D91", 
				"Wildtype" = "black")) + 
  xlab("Number of pauses (thousands)") + ylab("AUC") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(legend.position = "none") 

ggsave(plot = p1, filename = "AUCnPauseCor.svg", width = 5, height = 5)
