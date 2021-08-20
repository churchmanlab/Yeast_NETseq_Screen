#!/n/app/R/4.0.1/bin/Rscript

# Calculate % pauses in reads across deletion strains

# Load libraries
library(ggplot2)
library(svglite)
library(reshape2)

# Get arguments
args <- commandArgs(trailingOnly = TRUE)
t <- args[1]


# Load data
df = read.table(paste0("./pauseStrength_cov",t,".txt"), stringsAsFactors = FALSE)
colnames(df) = c("Mut", "nRead1", "nRead2", "nPause1", "nPause2")

# Calculate percent
df$Perc1 = (df$nPause1 / df$nRead1) * 100
df$Perc2 = (df$nPause2 / df$nRead2) * 100

# Calculate median
df$mPerc = apply(df[,c('Perc1', 'Perc2')], 1, median)

# Reorder for plotting
df$Mut = factor(df$Mut, levels = df[order(df$mPerc), 'Mut'])

# Plot
p = ggplot(df, aes(Mut, mPerc)) + 
  geom_bar(stat = "identity", fill = "forestgreen") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave(filename = paste0("pauseStrength_cov",t,".pdf"), device="pdf",plot = p, height = 6, width = 10)

# Correlate with sequencing depth
readMed = read.table("readTotalsMedians.txt", stringsAsFactors = FALSE)
colnames(readMed) = c("Mut", "ReadTotal")
readMed <- readMed[(readMed$Mut %in% df$Mut),] # Remove strains below gene number threshold

# Merge
data = data.frame(cbind(as.vector(df$Mut), df$mPerc, readMed$ReadTotal))
colnames(data) = c("Mut", "mPerc", "ReadTotal")

c = cor(df$mPerc, readMed$ReadTotal, method  = 'pearson')
cat(paste("There is a ", round(c^2, 4), " pearson correlation between pause strength and median sequencing depth.\n", sep=""))

# Plot
p = ggplot(data, aes(x = as.numeric(ReadTotal)/1000000, y = as.numeric(mPerc), label = Mut)) + 
  geom_point() + 
  geom_text() + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave(filename = paste0("pauseStrengthCovCor_cov",t,".pdf"),device="pdf", plot = p, height = 5, width = 5)
