# DNA Methyl-seq Differential Methylation Analysis

## Background 

This project aimed to identify genes that were differentially methylated in the amygdala and hippocampus between two closely related inbred mouse strains, and between social strata of an idiopathic mouse model for autism spectrum disorder. DNA was extracted from tissue punches of the respective regions and libraries were generated using the Agilent SureSelect XT target-enrichment protocol and sequenced on a NovaSeq S1 full flowcell. 

#### Overall workflow
- QC-filtering/adapter trimming (cutadapt, FASTQC, bash)
- Alignment (Bismark - the bisulfite-converted primary assembly of GRCm38, bash)
- Methylation estimates (CpG, CHG, CHH) (Bismark Methylation Extractor, bash)
- Differential methylation analysis (RNbeads \ DMRseq, R)

## Description of contents
- [bismark_pipeline.sh](https://github.com/stephen-siecinski/social_behavior_in_mice/blob/main/methyl-seq/bismark_pipeline.sh)
  - Processes raw FASTQ files. Filters reads based on phred score, size, error rate, performs multi-QC to evaluate sequencing quality, aligns to the bisulfite-converted mouse reference genome (GRCm38). CpG, CpH, CHH methylation estimates generated using Bismark Methylation Extractor.
- [RNbeads-runconfig.Rmd](https://github.com/stephen-siecinski/social_behavior_in_mice/blob/main/methyl-seq/RNbeads-runconfig.Rmd)
  - Performs an array of post-alignment QC, qualitative visualizations, and position / region-specific analyses of differential methylation. 
- [RNbeads_interpretation_visualization.Rmd](https://github.com/stephen-siecinski/social_behavior_in_mice/blob/main/methyl-seq/RNbeads_interpretation_visualization.Rmd)
  - Additional visualizations of the output dataframes from the RNbeads pipeline
- [dmr-seq.Rmd](https://github.com/stephen-siecinski/social_behavior_in_mice/blob/main/methyl-seq/dmr-seq.Rmd)
  - An algorithm that uses regional smoothing and permutation-based statistical tests to identify differentially methylated regions between experimental groups. 
- [dmr_gff_annotation.rmd](https://github.com/stephen-siecinski/social_behavior_in_mice/blob/main/methyl-seq/dmr_gff_annotation.rmd)
  - Annotate the regions identified by DMRseq using gff files for regulatory features, motif features, and gene annotations. Summarize results. 

## Summary results

#### DMRseq output
A representative plot depicting a differentially methylated region in the hippocampus between C57BL/6J and C58/J mice. 

<img src="https://github.com/stephen-siecinski/social_behavior_in_mice/blob/main/methyl-seq/dmr.png" width="800">

#### RNbeads output
Euclidean clustering of hippocampus brain samples based on differentially methylated CpG sites.

<img src="https://github.com/stephen-siecinski/social_behavior_in_mice/blob/main/methyl-seq/mouse_methyl_cluster.png" width="800">
