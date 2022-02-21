#!/n/app/R/4.0.1/bin/Rscript

# Plot hierarchically-clustered heatmap for all GO terms enriched in DE genes in at least 1 deletion strain

# Written by: K. Lachance
# Date: May 8, 2019

# Use: 

# Load libraries
library(ggplot2)
library(svglite)
library(reshape2)
# library(data.table)

# Get data
args <- commandArgs(trailingOnly = TRUE)

set <- args[1] # UP DOWN ALL
mat <- read.table(paste0("GOoverlaps",set,".txt"), sep = "\t", stringsAsFactors = FALSE, row.names=NULL, header=FALSE)
mut_list <- read.table("deletionStrains.txt", stringsAsFactors = FALSE)$V1

num_go = nrow(mat)
num_mut = length(mut_list)

# Create a data frame that can hold GO x deletion strain information
df = data.frame(matrix(0, nrow = num_go, ncol = num_mut), stringsAsFactors = FALSE)
colnames(df) = mut_list
rownames(df) = mat[,1]

# Loop through each GO term
for (i in 1:(num_go)) {

    muts = unlist(strsplit(mat$V2[i]," "))
    
    # Put 1 if that GO term is enriched in that deletion strain (otherwise 0)
    for (match_mut in muts) {    	
	df[i, colnames(df) == match_mut] = 1
    } # End of inner for loop

} # End of outer for loop

# Write output for completeness
write.table(df, paste0("plotGO_",set,"_output.txt"), quote = F, sep = "\t", row.names = T, col.names = T)

# Melt for plotting
m_df = melt(as.matrix(df))

# Order with hierarchical clustering
ord <- hclust(dist(df), method="average")$order
m_df$Var1 <- factor(m_df$Var1, levels=rownames(df)[ord])
m_df$Var2 <- factor(m_df$Var2, levels=colnames(df)[ord])

# Plot
p1 <- ggplot(m_df, aes(x = Var2, y = Var1, fill = value)) + 
   geom_tile(size=.7) + 
   scale_fill_gradient2(name=NULL, low="white", high="darkorchid1", labels = NULL) +  # coral3 for down-regulated, dodgerblue3 for up-regulated genes
   theme_bw() + 
   theme(axis.title=element_blank(), axis.text.y=element_text(size=2), axis.ticks.y = element_blank()) +
   theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 4, hjust = 1))+
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
   theme(legend.position="none")
   
ggsave(file=paste0("GOheatmap",set,".pdf"), plot=p1, width=10, height=10)

