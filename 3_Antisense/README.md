# Antisense transcription
For analysis of divergent antisense transcription at tandem genes, the sum of reads on the strand opposite the coding gene from 100 bp upstream of the sense TSS to 600 bp upstream of the sense TSS is divided by the sum of reads from the sense TSS to 500 bp downstream of the TSS on the coding strand. The genes selected were tandem gene pairs with each gene transcribed in the same direction and were obtained from (Churchman and Weissman, 2011). For analysis of antisense transcription, Pol II genes that did not overlap with another coding gene were chosen. The region analyzed spanned from the TSS to the polyadenylation site as defined by taking the most abundant TSS and polyadenylation site from (Pelechano et al., 2013). The sum of reads from the antisense strand was divided by the sum of reads from the sense strand. For all analyses, the log2 antisense to sense ratio was used. To generate antisense heatmaps, the log2 RPKM of NET-seq reads was used. Analysis at coding genes ranged from 250 bp upstream of the TSS to 4000 bp downstream of the coding TSS. To allow comparison between mutant and wild-type samples, a pseudocount of 1 was added to every position in all samples before calculating the log2 RPKM. Differential heatmaps were calculated by taking the log2 ratio of mutant / wild-type RPKM at each position. The heatmap illustrating the correlation of antisense : sense transcription ratios across deletion strains was created by calculating Pearson correlation of this ratio across all protein-coding non-overlapping genes. However, to increase the illustrated dynamic range, the first principal component was computationally removed before calculating these correlations.

# To quantify antisense transcription in the gene body as well as upstream tandem to genes
1. 