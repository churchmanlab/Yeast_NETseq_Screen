# Extracting pause positions
Pauses were identified in previously annotated transcription units (Xu et al., 2009) of well-expressed genes (average of > 2 reads per base-pair in two replicates). Pauses were defined as having reads higher than 3 standard deviations above the mean of the surrounding 200 nucleotides which do not contain pauses. Mean and standard deviation were calculated from a negative binomial distribution fit to the region of interest. Pauses were required to have at least 4 reads regardless of the gene’s sequencing coverage. Pauses were considered reproducible and used in downstream analyses when the Irreproducible Discovery Rate (IDR) is > 1% between two replicates. To calculate the IDR of each pause, log10 of pause strength (number of reads in pause) for each replicate was used as a proxy for pause score. IDR was calculated using the est.IDR function of the idr R package (mu = 3, sigma = 1, rho = 0.9, p = 0.5) (Li et al., 2011) . 

# Pol II pausing loci
Pause density was calculated as the ratio of total number of pauses to the total length of the genome considered when extracting pause positions (combined length of all well-expressed genes in both replicates of each deletion strain). To identify deletions that induced similar pausing patterns, 8,816 pauses were found to be shared by at least 8 deletion strains. Shared pauses were visualized with a heatmap, clustered on both axes using the eisenCluster correlation clustering method in the hybridHclust R package (Chipman and Tibshirani, 2006), which takes into account missing data (where there was not enough coverage to confidently identify pausing in a particular deletion strain). When calculating distribution of pauses across the gene body, all genes in which pauses were identified were normalized in length; the 5’ gene region was defined as the first 15% of each gene, the mid-gene region was defined as extending from the 15th percentile of gene length to the 85th percentile, and the 3’ gene region was defined as starting at 85% of gene length and extending to the annotated poly(A) site. The scrambled control for the pausing location analysis was created by randomly scrambling all identified pauses in all deletion strains across the gene in which they were discovered. All analyses related to sequence motifs underlying pause loci were conducted using the MEME suite of tools (Bailey et al., 1994, 2009).  The sequence ± 10 bp around each identified, reproducible pause (as well as the matched scrambled control) was extracted and used to run the MEME tool using parameters to find 0 - 1 motif per sequence, motifs 6 - 21 bp in length, and up to 10 motifs with an E-value significance threshold of 0.05 (Bailey et al., 1994). These significant motifs were compared to known transcription factor binding site motifs in the YEASTRACT_20130918 database (Teixeira et al., 2014) using the TOMTOM tool (Gupta et al., 2007) using default parameters, calling all hits as significant with an E-value greater than 0.1. TOMTOM searches were only performed on those motifs with a relative entropy greater than 5 and only the top match is reported.

# To extract pause loci

Before proceeding with this section, the "stt" bedGraph files created in 1_AlignData need to be sorted by chromosome. To do this, go back to that directory and run 
```
reps="yfg1-1 yfg1-2 goi1-1 goi1-2 wt-1 wt-2"
for rep in reps
do
sbatch ./Scripts/SortBedGraph.sh ${rep}_stt_pos.bedGraph ${rep}_stt_neg.bedGraph
done
```

1. First, pause loci need to be identified in each replicate. This can be done using 2 scripts consecutively: first find genes that pass the coverage threshold (which is input in the command line) with the following code block (note sample names need to be changed):  

```
samples="yfg1 goi1 wt"
mkdir logs
t=2 # Coverage threshold (average reads per nt)
strands="pos neg" 
for strand in $strands
do
for sample in $samples
do
sbatch -o logs/5_findGenes_cov${t}_${sample}_${strand}.log -e logs/5_findGenes_cov${t}_${sample}_${strand}.err ./Scripts/run_findGenesOverCovThreshold.sh ${sample} $t $strand
done
done
```
This will result in the output files `${mut}.covThreshold${t}.${strand}.bed`.  
Then identify pauses in those genes (note sample names need to be changed):

```
samples="yfg1 goi1 wt"
t=2
strands="pos neg" #"pos neg"
for strand in $strands
do
for sample in $samples
do
sbatch -o logs/5_findPauses_cov${t}_${sample}-1_${strand}.log -e logs/5_findPauses_cov${t}_${sample}-1_${strand}.err ./Scripts/run_findPausesRepCov.sh $sample ${sample}-1 $t $strand

sbatch -o logs/5_findPauses_cov${t}_${sample}-2_${strand}.log -e logs/5_findPauses_cov${t}_${sample}-2_${strand}.err ./Scripts/run_findPausesRepCov.sh $sample ${sample}-2 $t $strand
done
done
```
This will result in the output files `{sample}-1_cov${t}.pausesNB_T3.w100.${strand}.bed` and `{sample}-2_cov${t}.pausesNB_T3.w100.${strand}.bed`. 


2. Next, calculate the overall reproducibility and identify those reproducible pauses using IDR (note sample names need to be changed):
```
echo -e 'yfg1\ngoi1\nwt' > deletionStrains.txt
t=2
sbatch -o logs/5_IDR_${t}.log -e logs/5_IDR_${t}.err ./Scripts/IDR.sh $t
```
This will produce a file of reproducible pauses, `SampleName_PASS_pauseScore_cov${t}.txt`,which can be transformed into a bedGraph file and used in a genome browser such as IGV to look at data similar to that displayed in **Figure 4B**. It will also produce a scatterplot, a file called `SampleName_IDRscatter_cov${t}.pdf`, which is similar to that shown in **Figure 4A**. Note that this requires a file called `readTotals.txt` with 2 columns, one for the name of each replicate, and the other with the total number of uniquely-mapping reads in that sample.  To make this, use `./Scripts/GetReadTotals_MC.sh`  

3. Reproducibility of pauses across all deletion strains can be assessed using the command `./Scripts/plot_IDRrep_MC.R $t`, which will produce the plot `IDRrep_cov${t}.pdf`, which is similar to **Figure S5A**  

4. To calculate pause density, as in **Figure 4C**, use
```
t=2
sbatch -o logs/5_calcPDbyGene_${t}.log -e logs/5_calcPDbyGene_${t}.err ./Scripts/calcPauseDensitybyGene.sh deletionStrainsCovFiltered.txt $t
```
This will produce files containing the pause density of every gene in every deletion strain (files called `SampleName_cov${t}.pauseDensity.txt`) as well as a box plot called `pauseDensity_cov${t}.pdf`  

5. To calculate pause strength (median percent of reads in pauses), as shown in **Figure S5B**, use 
```
t=2
sbatch -o logs/5_calcPauseStrength.log -e logs/5_calcPauseStrength.err ./Scripts/calcPauseStrength.sh deletionStrains.txt $t
```
to generate the plot `pauseStrength_cov${t}.pdf`.These values can be correlated to sequencing depth, as in **Figure S5C** with the same script, which also produces the plot `pauseStrengthCovCor_cov${t}.pdf`.

For the next two steps, all pause loci for all deletion strains must be converted into vector form. This can be done with
```
t=2
sbatch -o logs/5_makeMutVector.log -e logs/5_makeMutVector.err ./Scripts/run_makeMutVector_MC.sh deletionStrains.txt ${t}
```

6. In order to perform the PCA analysis on shared pause loci across deletion strains, as depicted in **Figure S5D**, run `./Scripts/pausePCA.R` to produce the final plot `pausePCA.pdf`

7. Shared pause loci can also be visualized as a heatmap, like that shown in **Figure 4D**. In order to generate this figure, 
run `./Scripts/plot_pauseLociHeatmap.R`, which will generate the plot saved to a file called `pauseHeatmap.pdf`. Update line 22 in this script before running to require more or fewer strains to share a pause (currently set at 8). Note that this script may take a while to run!  


**Following steps are not yet updated**  
```
t=2
sbatch -o logs/5_makeLociVector.log -e logs/5_makeLociVector.err ./Scripts/makeLociVector.sh deletionStrains.txt ${t}
```

8. Pol II pause position varies across the gene body, as illustrated in **Figure 4E**. Before this figure can be reproduced, first all pauses must be shuffled within high-coverage genes. This can be done with the command `./Scripts/shufflePauses.sh SampleName`, which will prouce the file `SampleName.IDRrep.SHUFFLE.bed`. Combine all shuffled pause files together with the command `cat *.IDRrep.SHUFFLE.bed > ALL.IDRrep.SHUFFLE.bed`.  

9. Now, real and shuffled pauses can be overlapped with regions of the gene body, which can be done with the command `./Scripts/countOverlap.sh deletionStrains.txt`. Finally, we can generate the bar plot with proportion of pauses in each gene region for each deletion strain with the command `./Scripts/plotOverlap.R`, which will produce the file `PauseGeneOverlap.svg`
