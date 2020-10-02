#!/n/app/R/3.6.1/bin/Rscript

# Written by: K. Lachance
# Date: July 6, 2020

# Box plot of pausing indices

# Use: ./plot_boxPlot.R

# Load libraries
library(ggplot2)
library(svglite)

# Read in tables
data <- read.table("ALL.PI.bed", stringsAsFactors = FALSE)
colnames(data) = c("Chr", "Start1", "End1", "Start2", "End2", "Gene", "Skip", "Strand", "ROIcount", "notROIcount", "ROInorm", "notROInorm", "PI", "Region", "Mut")

data$Region = factor(data$Region, levels = c("TSS", "pA", "5pSS", "3pSS", "Anti"))

m = aggregate(data[data$Region == "TSS", 'PI'], by=list(Mut = data[data$Region == "TSS", 'Mut']), FUN=function(x) median(x, na.rm = TRUE))

write.table(m, "TSS_medians.txt", sep="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

ord = m[order(m$x),'Mut']
data$Mut = factor(data$Mut, levels = ord)

data$COL = "A"
data[data$Mut == "wt", 'COL'] = "B"

# Add coloring label
wt_median = m[(m$Mut == 'wt'), 'x']

wt_lower = quantile(data[(data$Mut == 'wt' & data$Region == "TSS"), 'PI'], 0.45, na.rm = TRUE)
wt_upper = quantile(data[(data$Mut == 'wt' & data$Region == "TSS"), 'PI'], 0.55, na.rm = TRUE)

# Plot
p1 <- ggplot(data[data$Region == "TSS",], aes(x = Mut, y = PI, color = COL, fill = COL)) + 
  geom_boxplot(outlier.shape = NA, coef = 0, notch = TRUE) + 
  geom_hline(aes(yintercept=wt_median), color="black") +
  geom_hline(aes(yintercept=wt_lower), color="black", linetype = "dotted") +
  geom_hline(aes(yintercept=wt_upper), color="black", linetype = "dotted") +
  coord_cartesian(ylim = c(0, 5)) + 
  scale_color_manual(values = c("A" = "#006838", "B" = "black")) + 
  scale_fill_manual(values = c("A" = "white", "B" = "#006838")) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(legend.position = "none") + 
  theme(strip.background = element_blank(),strip.text = element_blank(), panel.background = element_blank())

ggsave(filename = "PI_TSS_boxPlot.svg", plot = p1, height = 6, width = 10)

m = aggregate(data[data$Region == "pA", 'PI'], by=list(Mut = data[data$Region == "pA", 'Mut']), FUN=function(x) median(x, na.rm = TRUE))
write.table(m, "PA_medians.txt", sep="\t", quote = FALSE, row.names = FALSE, col.names =  FALSE)
ord = m[order(m$x),'Mut']
data$Mut = factor(data$Mut, levels = ord)

# Add coloring label
wt_median = m[(m$Mut == 'wt'), 'x']

wt_lower = quantile(data[(data$Mut == 'wt' & data$Region == "pA"), 'PI'], 0.45, na.rm = TRUE)
wt_upper = quantile(data[(data$Mut == 'wt' & data$Region == "pA"), 'PI'], 0.55, na.rm = TRUE)

p2 <- ggplot(data[data$Region == "pA",], aes(x = Mut, y = PI, color = COL, fill = COL)) +
  geom_boxplot(outlier.shape = NA, coef = 0, notch = TRUE) +
  geom_hline(aes(yintercept=wt_median), color="black") +
  geom_hline(aes(yintercept=wt_lower), color="black", linetype = "dotted") +
  geom_hline(aes(yintercept=wt_upper), color="black", linetype = "dotted") +
  coord_cartesian(ylim = c(0, 1.5)) +
  scale_color_manual(values = c("A" = "#BE1E2D", "B" = "black")) +
  scale_fill_manual(values = c("A" = "white", "B" = "#BE1E2D")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  theme(strip.background = element_blank(),strip.text = element_blank(), panel.background = element_blank())

ggsave(filename = "PI_pA_boxPlot.svg", plot = p2, height = 6, width = 10)


m = aggregate(data[data$Region == "5pSS", 'PI'], by=list(Mut = data[data$Region == "5pSS", 'Mut']), FUN=function(x) median(x, na.rm = TRUE))
write.table(m, "FiveSS_medians.txt", sep="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
ord = m[order(m$x),'Mut']
data$Mut = factor(data$Mut, levels = ord)

# Add coloring label
wt_median = m[(m$Mut == 'wt'), 'x']
wt_lower = quantile(data[(data$Mut == 'wt' & data$Region == "5pSS"), 'PI'], 0.45, na.rm = TRUE)
wt_upper = quantile(data[(data$Mut == 'wt' & data$Region == "5pSS"), 'PI'], 0.55, na.rm = TRUE)

p3 <- ggplot(data[data$Region == "5pSS",], aes(x = Mut, y = PI, color = COL, fill = COL)) +
  geom_boxplot(outlier.shape = NA, coef = 0, notch = TRUE) +
  geom_hline(aes(yintercept=wt_median), color="black") +
  geom_hline(aes(yintercept=wt_lower), color="black", linetype = "dotted") +
  geom_hline(aes(yintercept=wt_upper), color="black", linetype = "dotted") +
  coord_cartesian(ylim = c(0, 5)) +
  scale_color_manual(values = c("A" = "#2E3192", "B" = "black")) +
  scale_fill_manual(values = c("A" = "white", "B" = "#2E3192")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  theme(strip.background = element_blank(),strip.text = element_blank(), panel.background = element_blank())

ggsave(filename = "PI_5pSS_boxPlot.svg", plot = p3, height = 6, width = 10)

m = aggregate(data[data$Region == "3pSS", 'PI'], by=list(Mut = data[data$Region == "3pSS", 'Mut']), FUN=function(x) median(x, na.rm = TRUE))
write.table(m, "ThreeSS_medians.txt", sep="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
ord = m[order(m$x),'Mut']
data$Mut = factor(data$Mut, levels = ord)

# Add coloring label
wt_median = m[(m$Mut == 'wt'), 'x']
wt_lower = quantile(data[(data$Mut == 'wt' & data$Region == "3pSS"), 'PI'], 0.45, na.rm = TRUE)
wt_upper = quantile(data[(data$Mut == 'wt' & data$Region == "3pSS"), 'PI'], 0.55, na.rm = TRUE)

p4 <- ggplot(data[data$Region == "3pSS",], aes(x = Mut, y = PI, color = COL, fill = COL)) +
  geom_boxplot(outlier.shape = NA, coef = 0, notch = TRUE) +
  geom_hline(aes(yintercept=wt_median), color="black") +
  geom_hline(aes(yintercept=wt_lower), color="black", linetype = "dotted") +
  geom_hline(aes(yintercept=wt_upper), color="black", linetype = "dotted") +
  coord_cartesian(ylim = c(0, 5)) +
  scale_color_manual(values = c("A" = "#00AEEF", "B" = "black")) +
  scale_fill_manual(values = c("A" = "white", "B" = "#00AEEF")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  theme(strip.background = element_blank(),strip.text = element_blank(), panel.background = element_blank())

ggsave(filename = "PI_3pSS_boxPlot.svg", plot = p4, height = 6, width = 10)

m = aggregate(data[data$Region == "Anti", 'PI'], by=list(Mut = data[data$Region == "Anti", 'Mut']), FUN=function(x) median(x, na.rm = TRUE))
write.table(m, "Antisense_medians.txt", sep="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
ord = m[order(m$x),'Mut']
data$Mut = factor(data$Mut, levels = ord)

# Add coloring label
wt_median = m[(m$Mut == 'wt'), 'x']
wt_lower = quantile(data[(data$Mut == 'wt' & data$Region == "Anti"), 'PI'], 0.45, na.rm = TRUE)
wt_upper = quantile(data[(data$Mut == 'wt' & data$Region == "Anti"), 'PI'], 0.55, na.rm = TRUE)

p5 <- ggplot(data[data$Region == "Anti",], aes(x = Mut, y = PI, color = COL, fill = COL)) +
  geom_boxplot(outlier.shape = NA, coef = 0, notch = TRUE) +
  geom_hline(aes(yintercept=wt_median), color="black") +
  geom_hline(aes(yintercept=wt_lower), color="black", linetype = "dotted") +
  geom_hline(aes(yintercept=wt_upper), color="black", linetype = "dotted") +
  coord_cartesian(ylim = c(0, 5)) +
  scale_color_manual(values = c("A" = "#507BA7", "B" = "black")) +
  scale_fill_manual(values = c("A" = "white", "B" = "#507BA7")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  theme(strip.background = element_blank(),strip.text = element_blank(), panel.background = element_blank())

ggsave(filename = "PI_Antisense_boxPlot.svg", plot = p4, height = 6, width = 10)
