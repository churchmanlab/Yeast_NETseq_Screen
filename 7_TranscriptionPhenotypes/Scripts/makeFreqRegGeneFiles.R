#!/n/app/R/3.6.1/bin/Rscript

# Plot number of strains gene is DE in v. transcription phenotypes

# Written by: K. Lachance
# Date: August 28, 2020

# Use: ./makeFreqRegGenes.R UP/DOWN

library(ggplot2)
library(svglite)

args <- commandArgs(trailingOnly = TRUE)
DIR <- args[1]

# Get data
geneList = read.table(paste("DEgenes.", DIR, ".8orMore.txt", sep=""), stringsAsFactors = FALSE)$V1

ASdata = read.table("../3_Antisense/CodingAntisense/ALL.geneNETseqRatio.bed", stringsAsFactors = FALSE)
colnames(ASdata) <- c("Chr", "Start", "End", "Gene", "Skip", "Strand", "Sense", "Anti", "Ratio", "Mut")
wt_ASdata = ASdata[ASdata$Mut == "wt",]

TSSdata = read.table("../4_PausingIndex/ALL_TSS.PI.bed", stringsAsFactors = FALSE)
colnames(TSSdata) = c("Chr", "Start1", "End1", "Start2", "End2", "Gene", "Skip", "Strand", "ROIcount", "notROIcount", "ROInorm", "notROInorm", "PI", "Region", "Mut")
wt_TSSdata = TSSdata[TSSdata$Mut == "wt",]

# Create wildtype dataframes
wt_ass = wt_ASdata[wt_ASdata$Gene %in% geneList, c("Gene", "Ratio")]
colnames(wt_ass) = c("Gene", "WTvalue")
wt_tss = wt_TSSdata[wt_TSSdata$Gene %in% geneList, c("Gene", "PI")]
colnames(wt_tss) = c("Gene", "WTvalue")

# Make out DF
out_ass = data.frame("Gene" = character(), "WTvalue" = numeric(), "MUTvalue" = numeric(), "Mut" = character(), "Cat" = character(), stringsAsFactors = FALSE)
out_tss = data.frame("Gene" = character(), "WTvalue" = numeric(), "MUTvalue" = numeric(), "Mut" = character(), "Cat" = character(), stringsAsFactors = FALSE)

# Loop through each deletion strain and each gene
muts = read.table("deletionStrains.txt", stringsAsFactors = FALSE)$V1
for (mut in muts) {
    
    cat(paste("Working on ", mut, "!\n", sep=""))
    
    # Get data for mut
    mut_ass = ASdata[ASdata$Mut == mut,]
    mut_tss = TSSdata[TSSdata$Mut == mut,]

    # Make temporary DFs to hold info
    tmp_ass = merge(wt_ass, mut_ass, by = "Gene")[,c("Gene", "WTvalue", "Ratio")]
    tmp_ass$Mut = mut
    tmp_ass$Cat = "NO"
    tmp_tss = merge(wt_tss, mut_tss, by = "Gene")[,c("Gene", "WTvalue", "PI")]
    tmp_tss$Mut = mut
    tmp_tss$Cat = "NO"

    # Check if gene is DE in this strain
    DE_genes = read.table(paste("../2_DiffExpression/", mut, ".DEgenes.txt", sep=""), header = TRUE, row.names = 1, stringsAsFactors = FALSE)
    
    for (gene in tmp_ass$Gene) {
    	if (gene %in% rownames(DE_genes)) {
	   if (DIR == "UP" & DE_genes[gene,'log2FoldChange'] >= 1) {  # Up-regulated
	      tmp_ass[tmp_ass$Gene == gene, 'Cat'] = "YES"
	   } 
	   if (DIR == "DOWN" & DE_genes[gene,'log2FoldChange'] <= -1) {	# Down-regulated
	      tmp_ass[tmp_ass$Gene == gene, 'Cat'] = "YES"
	   }
	}
    }   


    for (gene in tmp_tss$Gene) {
        if (gene %in% rownames(DE_genes)) {
           if (DIR == "UP" & DE_genes[gene,'log2FoldChange'] >= 1) {	# Up-regulated
              tmp_tss[tmp_tss$Gene == gene, 'Cat'] = "YES"
           }           
	   if (DIR == "DOWN" & DE_genes[gene,'log2FoldChange'] <= -1) { # Down-regulated
              tmp_tss[tmp_tss$Gene == gene, 'Cat'] = "YES"
           }
        }
    }

    # Add to lengthening dataframes
    out_ass = rbind(out_ass, tmp_ass)
    out_tss = rbind(out_tss, tmp_tss)

} # End of outer for loop

# Write tables
write.table(out_ass, paste("DE_AntisenseRatio.", DIR, ".txt", sep=""), row.names = FALSE, col.names = FALSE, quote = FALSE, sep="\t")
write.table(out_tss, paste("DE_TssPI.", DIR, ".txt", sep=""), row.names = FALSE, col.names = FALSE, quote = FALSE, sep="\t")
