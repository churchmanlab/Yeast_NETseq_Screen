#!/n/app/R/3.6.1/bin/Rscript

# Plot correlation between number of DE genes with doubling time and sequencing depth

# Written by: K. Lachance
# Date: May 8, 2019

# Use: ./plotDEgeneCor.R

# Load libraries
library(ggplot2)
library(svglite)

# Read table
dat <- read.table("Screen_DEcor.txt", header = FALSE, stringsAsFactors = FALSE)
colnames(dat) <- c("Mut", "nReads", "DoublingTime", "nDE", "Category")

# Calculate correlations
c1 = cor(dat$DoublingTime, dat$nDE, method = "pearson")
cat(paste("There is a Pearson correlation of ", round(c1, 2), " between doubling time and the number of identified DE genes.\n", sep=""))

c2 = cor(dat$nReads, dat$nDE, method = "pearson")
cat(paste("There is a Pearson correlation of ", round(c2, 2), " between median sequencing depth and the	number of identified DE genes.\n", sep=""))

# Plot
p1 <- ggplot(dat, aes(x = DoublingTime, y = nDE, color = Category, label = Mut)) + 
  geom_point() + 
  geom_text() + 
  scale_color_manual(values = c("RNAprocessing" = "#39B54A",
  			        "HistoneVariant" = "#939598",
				"ChromatinModifer" = "#F15A29", 
				"ChromatinRemodeler" = "#1C75BC", 
				"TranscriptionElongation" = "#662D91", 
				"Wildtype" = "black")) + 
  xlab("Doubling time (min)") + ylab("# DE genes") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(legend.position = "none") 

ggsave(plot = p1, filename = "DEgeneCor_DoublingTime.svg", width = 5, height = 5)

p2 <- ggplot(dat, aes(x = nReads / 1000000, y = nDE, color = Category, label = Mut)) +
  geom_point() +
  geom_text() +
  scale_color_manual(values = c("RNAprocessing" = "#39B54A",
                                "HistoneVariant" = "#939598",
                                "ChromatinModifer" = "#F15A29",
                                "ChromatinRemodeler" = "#1C75BC",
                                "TranscriptionElongation" = "#662D91",
                                "Wildtype" = "black")) +
  xlab("Median sequencing depth (millions)") + ylab("# DE genes") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none")

ggsave(plot = p2, filename = "DEgeneCor_SeqDepth.svg", width = 5, height = 5)
