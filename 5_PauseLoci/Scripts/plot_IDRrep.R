#!/n/app/R/3.6.1/bin/Rscript

# Plot reproducibility of p=3 pauses assessed with IDR

# Written by: K. Lachance
# Date: Feb 27, 2020

# Use: ./plot_IDRrep.R

# Load libraries
library(ggplot2)
library(reshape2)
library(svglite)

# Read data
data = read.table("IDR_reproducibility.txt", stringsAsFactors = FALSE)
colnames(data) = c("Mut", "Rep1", "Rep2", "Comb", "IDR", "Perc")

# Calculate non-reproducible pause number
data$NoComb = data$Rep1 + data$Rep2 - (2 * data$Comb)
data$noIDR = data$Comb - data$IDR

# Make data frame for plotting
df_idr = cbind(data$Mut, data$IDR, "zzzIDR")
df_comb = cbind(data$Mut, data$noIDR, "zComb")
df_nocomb = cbind(data$Mut, data$NoComb, "NoComb")

df = data.frame(rbind(df_idr, df_comb, df_nocomb), stringsAsFactors = FALSE)
colnames(df) = c("Mut", "N", "Lab")
df$N = as.numeric(df$N)

# Order by Comb for plotting
df$Mut <- factor(df$Mut, levels = df$Mut[order(data$IDR)])

# Plot

p <- ggplot(df, aes(x = Mut, y = N/1000, fill = Lab)) + 
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(values=c("zzzIDR" = "darkcyan", "zComb" = "cyan3", "NoComb" = "gray70")) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "none")

ggsave(file="IDRrep.svg", plot=p, width=10, height=6)

