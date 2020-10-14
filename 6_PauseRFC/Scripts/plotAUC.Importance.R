#!/n/app/R/3.6.1/bin/Rscript

# Plot heatmaps with importance and AUC values

# Written by: K. Lachance
# Date: May 11, 2020

# Use: ./plotAUC.Importance.R

# Import libraries
library(ggplot2)
library(reshape2)

# Import importance values
dat = read.table("ALL.importance.txt", stringsAsFactors = FALSE)

feature_list = c("PauseNucA", "PauseNucE", "PauseNucAB", "PauseNucBC", "PauseNucCD", "PauseNucDE", "pCG", "ShapeHelT", "ShapeMGW", "ShapeProT", "ShapeRoll", "GenePerc", "DistFromTSS", "h3k36me3", "h3k79me", "h3k79me3", "h4r3me2s", "h4k20me", "TmNN")
mut_list = unique(dat$V3)
df = data.frame(matrix(0, nrow = length(feature_list), ncol = length(mut_list)), stringsAsFactors = FALSE)
row.names(df) = feature_list
colnames(df) = mut_list

for (f in feature_list) {
  for (m in mut_list) {
    df[f, m] = dat[dat$V1 == f & dat$V3 == m, 2]
  }
}

dat = df

col.order <- hclust(dist(t(dat)))$order
COL.ORDER = col.order

dat_new <- dat[, col.order] # re-order matrix accoring to clustering

df_molten_dat <- melt(as.matrix(dat_new)) # reshape into dataframe
names(df_molten_dat)[c(1:2)] <- c("Feature", "Mut")

p = ggplot(data = df_molten_dat,aes(x = Mut, y = Feature, fill = value)) + 
  geom_raster() +
  scale_y_discrete(limits = rev(levels(df_molten_dat$Feature))) + 
  scale_fill_gradient2(low="black", high= "#FFC000", midpoint = 50) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text.y = element_blank()) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

ggsave(filename = "ImportanceHeatmap.svg", plot = p, width = 10, height = 6)

# Do the same for AUC values
dat = read.table("AUCvalues.txt",stringsAsFactors = FALSE)
colnames(dat) = c("Feature", "Mut", "value")
dat$Mut = factor(dat$Mut, levels = mut_list[COL.ORDER])

col.order <- hclust(dist(t(dat)))$order
dat_new <- dat[, col.order] # re-order matrix accoring to clustering

df_molten_dat <- melt(as.matrix(dat_new)) # reshape into dataframe
names(df_molten_dat)[c(1:2)] <- c("Feature", "Mut")

feature_list = c("Seq", "Shape", "Gene", "ChIP", "CTD", "COMB")
df_molten_dat$Feature = factor(df_molten_dat$Feature, levels = feature_list)

p = ggplot(data = df_molton_dat, aes(x = Mut, y = Feature, fill = value, label = round(value, 2))) + 
  geom_raster() +
  geom_text() + 
  scale_fill_gradient(low = "white", high = "darkorange") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

ggsave(filename = "AUC_heatmap.svg", plot = p, width = 6, height = 4)
