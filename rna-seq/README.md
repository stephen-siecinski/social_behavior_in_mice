# RNA-seq Differential Expression Analysis

## Background 

This project aimed to identify genes that were differentially expressed in the amygdala and hippocampus between two closely related inbred mouse strains, and between social strata of an idiopathic mouse model for autism spectrum disorder. RNA was extracted from tissue punches of the respective regions and libraries were generated using the Illumina TrueSeq mRNA LP kit and sequenced on a NovaSeq 6000 S1 full flowcell. 

#### Overall workflow
- QC-filtering/adapter trimming (cutadapt, FASTQC, bash)
- Alignment (STAR - GRCm38.96, bash)
- Generate counts for known transcripts (FeatureCount, bash)
- Principal component analysis (PCAtools, R)
- Differential expression analysis (DESeq2, R)
- Gene set enrichment analysis (fGSEA, R)

## Description of contents
- C58_C57_expanded_models_results.Rmd
- fGSEA_rnaseq.Rmd
- DESeq2_heatmaps.Rmd

