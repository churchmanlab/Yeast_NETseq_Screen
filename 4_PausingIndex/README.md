# Pausing index calculation
Pausing indices were calculated as the length-normalized Pol II density in the region of interest (-50 bp to +150 bp around TSS, ±100 bp around poly(A) sites, and ±10 bp around 5’ and 3’ splice sites) divided by the length-normalized Pol II density in the remainder of the gene. 

# To calculate pausing indices around RNA processing sites
1. Calculate pausing indices for each deletion strain with the command `./Scripts/run_calcPI.R`, which produces the following 5 files:
   - `SampleName_TSS.PI.bed`
   - `SampleName_pA.PI.bed`
   - `SampleName_Antisense.PI.bed`
   - `SampleName_3pSS.PI.bed`
   - `SampleName_5pSS.PI.bed`
2. To calculate the correlation between PIs within a single deletion strain, as shown in **Figure S3A**, run the command `./Scripts/plotGeneCor.R SampleName`, which will generate plots named `SampleName_...v..._PIscatter.svg`
3. Create heatmaps for each gene in each deletion strain, as in **Figure 3B-C**. This can be done by running the script `./Scripts/plot_geneHeatmap_geneEnds.R SampleName` to produce plots `SampleName_TSS_GenePIheatmap.svg`, `SampleName_pA_GenePIheatmap.svg`, and `SampleName_Antisense_GenePIheatmap.svg`. To do the same for splice sites, run the script `./Scripts/plot_geneHeatmap_splice.R SampleName` to produce plots `SampleName_5pSS_GenePIheatmap.svg` and `SampleName_3pSS_GenePIheatmap.svg`
4. To illustrate the cumulative density for all regions for a single sample, as in **Figure S3C**, use the command `./Scripts/plotCDFpi.R` to generate the plot `PI_CDF.svg`
5. To create heatmaps for PI similar to those in **Figure S3D-H**, you first need to generate a file containing all PI values across all deletion strains with the command `cat *.PI.bed > ALL.PI.bed`. Then, run the script `./Scripts/plot_geneHeatmap_all.R TSS/pA/Anti/5pSS/3pSS` to geneate plots named `ALL_..._GenePIheatmap.svg`
6. In order to plot the boxplot of all PI distributions across all deletion strains as in **Figure 3F-J**, use the command `./Scripts/plot_boxPlot.R`, which will produce plots called `PI_..._boxPlot.svg` as well as files with all of the median values for each region PI
7. To calculate the correlation between median PI values across all deletion strains, as in **Figure 3D-E** and **Figure S3B**, use the command `./Scripts/plotMedianCor.R` to generate plots with the names `Medians_...v..._PIscatter.svg`
8. Metagene plots, as shown in **Figure 3K-Q** and **Figure S4A-D**, require some pre-script work to generate the necessary files for plotting. First, run the command `./Scripts/run_overlap_data.sh`, which generates the overlap files . From here, you can run `./Scripts/run_metagene_plot.sh`, which generates files ready for plotting in the directory `readyToPlot`. Finally, you can run the script `./Scripts/metagenePlot.R SampleName` to generate all metagene plots (TSS, pA, antisense, and both splice sites) for a given sample, which will be saved to files called `SampleNamevWT_metagene_....svg`

# To calculate splicing indecies
1. 


