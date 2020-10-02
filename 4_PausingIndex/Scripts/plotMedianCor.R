#!/n/app/R/3.6.1/bin/Rscript

# Written by: K. Lachance
# Date: July 20, 2020

# Correlation between PIs in different regions across deletion strains

# Use: ./plotMedianCor.R

# Load libraries
library(ggplot2)
library(svglite)

# Read in tables
COLORS <- read.table("Mut_colors.txt", stringsAsFactors = FALSE)
colnames(COLORS) = c("Mut", "x")

data_TSS <- read.table("TSS_medians.txt", stringsAsFactors = FALSE)
colnames(data_TSS) = c("Mut", "Med")

data_pA <- read.table("PA_medians.txt", stringsAsFactors = FALSE)
colnames(data_pA) = c("Mut", "Med")

data_Anti <- read.table("Antisense_medians.txt", stringsAsFactors = FALSE)
colnames(data_Anti) = c("Mut", "Med")

data_5pSS <- read.table("FiveSS_medians.txt", stringsAsFactors = FALSE)
colnames(data_5pSS) = c("Mut", "Med")

data_3pSS <- read.table("ThreeSS_medians.txt", stringsAsFactors = FALSE)
colnames(data_3pSS) = c("Mut", "Med")

# Merge
start_end = merge(data_TSS, data_pA, by = "Mut")[,c("Mut", "Med.x", "Med.y")]
start_anti = merge(data_TSS, data_Anti, by = "Mut")[,c("Mut", "Med.x", "Med.y")]
end_anti = merge(data_pA, data_Anti, by = "Mut")[,c("Mut", "Med.x", "Med.y")]
splice_splice = merge(data_5pSS, data_3pSS, by = "Mut")[,c("Mut", "Med.x", "Med.y")]

# Correlation
c1 = cor(start_end$Med.x, start_end$Med.y, method = "pearson", use = "pairwise.complete.obs")
cat(paste("There is a Pearson correlation of ", round(c1, 2), " between TSS PI and pA PI medians across deletion strains.\n", sep=""))

c2 = cor(start_anti$Med.x, start_anti$Med.y, method = "pearson", use = "pairwise.complete.obs")
cat(paste("There is a Pearson correlation of ", round(c2, 2), " between TSS PI and antisense PI medians across deletion strains.\n", sep=""))

c3 = cor(end_anti$Med.x, end_anti$Med.y, method = "pearson", use = "pairwise.complete.obs")
cat(paste("There is a Pearson correlation of ", round(c3, 2), " between pA PI and antisense PI medians across deletion strains.\n", sep=""))

c4 = cor(splice_splice$Med.x, splice_splice$Med.y, method = "pearson", use = "pairwise.complete.obs")
cat(paste("There is a Pearson correlation of ", round(c4, 2), " between 5pSS PI and 3pSS PI medians across deletion strains.\n", sep=""))

# Merge for color
start_end2 = merge(start_end, COLORS, by = "Mut")
start_anti2 = merge(start_anti, COLORS, by = "Mut")
end_anti2 = merge(end_anti, COLORS, by = "Mut")
splice_splice2 = merge(splice_splice, COLORS, by = "Mut")

# Plot
p1 <- ggplot(start_end2, aes(x = Med.x, y = Med.y, color = x, label = Mut)) +
  geom_point() +
  geom_text() + 
  scale_color_manual(values = c("A" = "#39B54A", "B" = "#939598", "C" = "#F15A29", "D" = "#1C75BC", "E" = "#662D91", "F" = "black")) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())	+ 
  theme(legend.position = "none")

ggsave(filename = "Medians_TSSvPA_PIscatter.svg", plot = p1, height = 4, width = 4)

p2 <- ggplot(start_anti2, aes(x = Med.x, y = Med.y, color = x, label = Mut)) +
  geom_point() +
  geom_text() +
  scale_color_manual(values = c("A" = "#39B54A", "B" = "#939598", "C" = "#F15A29", "D" = "#1C75BC", "E" = "#662D91", "F" = "black")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none")

ggsave(filename = "Medians_TSSvAnti_PIscatter.svg", plot = p2, height = 4, width = 4)

p3 <- ggplot(end_anti2, aes(x = Med.x, y = Med.y, color = x, label = Mut)) +
  geom_point() +
  geom_text() +
  scale_color_manual(values = c("A" = "#39B54A", "B" = "#939598", "C" = "#F15A29", "D" = "#1C75BC", "E" = "#662D91", "F" = "black")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none")

ggsave(filename = "Medians_pAvAnti_PIscatter.svg", plot = p3, height = 4, width = 4)

p4 <- ggplot(splice_splice2, aes(x = Med.x, y = Med.y, color = x, label = Mut)) +
  geom_point() +
  geom_text() +
  scale_color_manual(values = c("A" = "#39B54A", "B" = "#939598", "C" = "#F15A29", "D" = "#1C75BC", "E" = "#662D91", "F" = "black")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none")

ggsave(filename = "Medians_5pSSv3pSS_PIscatter.svg", plot = p4, height = 4, width = 4)
