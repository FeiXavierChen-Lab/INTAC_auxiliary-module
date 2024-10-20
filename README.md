# Catalytic-independent functions of INTAC confer sensitivity to BET inhibition
Primary sequencing data and BigWig files are deposited at the GEO depository: GSE274476.

**RNA-seq analysis**

The paired raw reads were quality-checked with FastQC v0.11.9 (Babraham Institute) and trimmed by Trim Galore v0.6.6 (Babraham Institute) with the parameter “-q 25” to remove adaptors and low-quality reads. The remaining reads were aligned to the mouse mm10 and human hg19 assemblies using STAR v2.7.5c (Dobin et al., 2013) with parameters “--outSAMtype BAM SortedByCoordinate –twopassModeBasic --outFilterMismatchNmax  3”. Duplicates and low-quality reads were removed, and reads mapped in proper pairs were isolated using SAMtools v1.9102 with parameter “-f 2”. The number of spike-in hg19 reads was counted with SAMtools v1.9 and used to generate normalization factor alpha = 1e6/ hg19_count for coverage profiles. Strand-specific normalized bigwigs were generated.  Raw gene counts were generated by featureCounts tool from the Rsubread R package v2.0.1 (Liao et al., 2019). Differential expression analysis was performed using DESeq2 R package v1.26.0 (Love et al., 2014), and the counts for the human hg19 spike-in were used to estimate the size factors. Genes with a false discovery rate < 0.05 and absolute log2FC > 1 were considered as differentially expressed genes.

**CUT&Tag analysis**

Raw CUT&Tag reads were trimmed using Trim Galore v.0.6.6 (Babraham Institute) in paired-end mode. Trimmed reads were aligned to mouse mm10 genome using Bowtie v.2.4.4 with the parameters “-N 1 -L 25 -X 700 --no-mixed --no-discordant” (Langmead and Salzberg, 2012). Duplicate reads were removed with Picard Tools v.2.25.5 (Broad Institute) and the reads were shifted to compensate for the offset in tagmentation site relative to the Tn5 binding site using the alignmentSieve function of deepTools v.3.5.1 with the ‘--ATACshift’ option (Ramirez et al., 2016). 
