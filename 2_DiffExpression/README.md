# Differential gene expression and gene ontology enrichment analysis
Expression levels were normalized using DESeq2 (Love et al., 2014). Differential expression analysis between deletion strains (with replicates) and wild-type strains was performed using DESeq2 (Love et al., 2014) for all genes in Xu et al 2009, as well as their corresponding "antisense genes" which were defined as the opposite strand for all genes, where there was not overlap with another known genes. Genes were considered differentially expressed if they had an adjusted p-value < 0.05 and an absolute log2 fold change > 1.0. GO term enrichment analysis was performed using Ontologizer (http://ontologizer.de) using input files downloaded from the GO website (http://geneontology.org): go.obo and sgd.gaf. The .gaf association file was edited to add a new GO term "antisense transcripts" associated with all antisense genes as defined above.  

# To identify differentially expressed genes
1. Create files required for running DESeq2 (counts and sample matrix) by running `./Scripts/makeDESeqMatrix.sh Sample-1 Sample-2 SampleName`. This script can accommodate up to 4 wildtype samples and 6 perturbation samples. Running this command will result in two files:
   - `SampleName.mat.txt`, which is a gene x sample matrix containing un-normalized gene expression calculated from the full NET-seq read
   - `SampleName.cond.txt`, which is a properly formatted list of samples
2. Run DESeq2 with the command `./Scripts/DESeqAnalysis.R SampleName`. This will report the total number of DE genes and produce another two files:
   - `SampleName.ALLgenes.txt`, which is the complete output from running DESeq2
   - `SampleName.DEgenes.txt`, which is the output of DESeq2, filtered to only include those genes with pAdj < 0.05 and |Log2FoldChange| > 1. This can be used to create tables similar to **Table S2**
3. Create a volcano plot illlustrating differentially expressed genes with the command `./Scripts/plotVolcano.R SampleName`, which will save the plot to a file called `SampleName_volcano.pdf`. Before running this script, optionally remove "anti" genes from SampleName.ALLgenes.txt: `for mut in $muts;
do
grep -v '^anti' ${mut}.ALLgenes.txt > ${mut}.ALLgenesNoAnti.txt;
done`
5. To re-create **Figure 1B**, run the command `./Scripts/plotNumDEgenes.R` to produce the plot called `NumDEgenes.pdf`. The input file for this script, `Screen_DE.txt` was created by transcribing the output of the volcano plot script (Step 3) for each deletion strain (4 columns: Deletion strain name, number of differentially up-regulated genes, number of differentially down-regulated genes, total number of differentially expressed genes): `touch Screen_DE.txt
for mut in $muts;
do
awk '{print $1"\t"$3"\t"$7"\t"$14}' logs/volcano_${mut}_noAnti.log >> Screen_DE.txt;
done`
6. To re-create **Figure 1C**, run the command `./Scripts/plotCDFforDEgenes.R` to produce the plot called `CDFforDEgenes.pdf`. The input file for this script, `DEgeneCount.txt`, was created by running the following command after DESeq2 (Step 2) was complete for each deletion strain (optionally first remove "anti" genes from SAmpleName.DEgene.txt as above): `cat *DEgenesNoAnti.txt | cut -f 1 | grep -v baseMean | sort | uniq -c | awk 'BEGIN {OFS="\t"} {print $1, $2}' | sort -k1 -n  > DEgeneCount.txt`. A similar process can be completed for the enriched GO terms (see below) to produce **Figure S1E**
7. To re-create **Figure S1C-D**, run the command `./Scripts/plotDEgeneCor.R` to produce the plots called `DEgeneCor_DoublingTime.svg` and `DEgeneCor_SeqDepth.svg`. The input file fot this script, `Screen_DEcor.txt`, was created such that each row contains 5 columns (Deletion strain name, median number of reads sequenced across replicates, median doubling time of deletion strain in minutes, total number of DE genes observed, and the category of the factor deleted, as in **Figure 1B**). 

All of these analyses may be run while only considering the interior gene body to avoid potential biases of Pol II pausing at the TSS or pA sites. This analysis (which produced **Figure S1B**) is identical to the above analysis, but Step 1 is run with the command `./Scripts/makeDESeqMatrix_subgenic.sh Sample-1 Sample-2 SampleName` instead. 

# To identify gene ontology terms enriched in differentially expressed genes
1. Create lists of differentially expressed genes (just gene name, split by whether that gene is up or down regulated, or keep both together). This can be done with the following commands:
`for mut in $muts;
do`  
`cat ${mut}.DEgenes.txt | grep -v baseMean | awk '($3 > 1) {print $0}' | cut -f1 > ${mut}.DEgeneList.UP.txt`  
`cat ${mut}.DEgenes.txt | grep -v baseMean | awk '($3 < 1) {print $0}' | cut -f1 > ${mut}.DEgeneList.DOWN.txt`    
`cat ${mut}.DEgenes.txt | grep -v baseMean | cut -f1 > ${mut}.DEgeneList.ALL.txt`  
`done`
2. Download ontologizer.jar (http://ontologizer.de/commandline/)
3. Perform gene ontology enrichment analysis by running  `for mut in $muts; do java -jar /pathToProgram/Ontologizer.jar -a ../0_Annotations/GO_databases/sgd_mtc.gaf -g ../0_Annotations/GO_databases/go.obo -s ${mut}.DEgeneList.UP/DOWN/ALL.txt -p ../0_Annotations/Genes/txpts_all_and_antisense_GENES.txt -c Parent-Child-Intersection; done`  This can be run in turn for each of the UP, DOWN, and ALL DEgeneLists.
4. Parse tables, filter and make gene lists with ParseGOtable.R: `for mut in $muts;
do
sbatch -p short -t 0-00:20 --wrap="../Scripts/ParseGOtable.R $mut Parent-Child-Intersection";
done`

5. Create a list of all GO terms that are enriched in at least one deletion strain, as `ALL.GO.txt`, e.g. `cat *_UP_GOterms.txt | sort -u > ALL.UP_GO.txt
`
6. Create a matrix for heatmap creation with the command `./Scripts/findCommonGO_MC.sh`. This will produce the file `GOoverlaps.txt`
7. Plot heatmap as in **Figure 1E** with the command `./Scripts/plotGOheatmap_MC.R UP/DOWN/ALL` to produce the plot called `GOheatmap.pdf` and the additional file `plotGO_output.txt`.  
8. Run MakeGOtables.R to produce tables similar to those in **Table S3**
