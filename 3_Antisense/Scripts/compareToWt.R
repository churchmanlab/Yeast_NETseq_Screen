#!/n/app/R/3.6.1/bin/Rscript

# Compare log2(Sense) and log2(Antisense) for each gene, mutant v WT

# Written by: K. Lachance
# Date: October 30, 2018

# Use: ./sense_antisense_compareToWt.R ALL.NETseqRatio.bed coding / tandem

# Read in parameters (data)
args <- commandArgs(trailingOnly = TRUE)
dat_file <- args[1]
label <- args[2]

# Read table
mat <- read.table(args[1], stringsAsFactors = FALSE)
colnames(mat) <- c("Chr", "Start", "End", "Gene", "Value", "Strand", "Sense", "Antisense", "Ratio", "Mutant")

# Loop through each gene / mut combo
for (i in 1:nrow(mat)) {
    
    if (i %% 100 == 0) {cat(paste("Done with ", i, " / ", nrow(mat), "rows!\n"))}

    g = mat$Gene[i]
    wt_sense = mat[mat$Gene == g & mat$Mutant == "wt", 'Sense']
    wt_anti = mat[mat$Gene == g & mat$Mutant == "wt", 'Antisense']

    mat$SenseWTratio[i] = mat$Sense[i] / wt_sense
    mat$AntiWTratio[i] = mat$Antisense[i] / wt_anti

}

# Remove wt
out = mat[mat$Mut != "wt", ]

# Get file name
if (label == "coding") {
   filename = "CodingAntisense/ALL.geneNETseqRatio_withWTratio.bed"
} else {
   filename = "TandemAntisense/ALL.tandemNETseqRatio_withWTratio.bed"
}

# Write table
write.table(out, "ALL.geneNETseqRatio_withWTratio.bed", quote = F, row.names = F, col.names = F, sep = "\t")