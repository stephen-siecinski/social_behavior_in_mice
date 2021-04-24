# DNA Methyl-seq Differential Methylation Analysis

## Background 

This project aimed to identify genes that were differentially methylated in the amygdala and hippocampus between two closely related inbred mouse strains, and between social strata of an idiopathic mouse model for autism spectrum disorder. DNA was extracted from tissue punches of the respective regions and libraries were generated using the Agilent SureSelect XT target-enrichment protocol and sequenced on a NovaSeq S1 full flowcell. 

#### Overall workflow
- QC-filtering/adapter trimming (cutadapt, FASTQC, bash)
- Alignment (Bismark - the bisulfite-converted primary assembly of GRCm38, bash)
- Methylation estimates (CpG, CHG, CHH) (Bismark Methylation Extractor, bash)
- Differential methylation analysis (RNbeads \ DMRseq, R)

## Description of contents
- [bismark_pipeline.sh](https://github.com/stephen-siecinski/social_behavior_in_mice/blob/main/methyl-seq/bismark_pipeline.sh):
- [RNbeads-runconfig.Rmd](https://github.com/stephen-siecinski/social_behavior_in_mice/blob/main/methyl-seq/RNbeads-runconfig.Rmd):
- [RNbeads_interpretation_visualization.Rmd](https://github.com/stephen-siecinski/social_behavior_in_mice/blob/main/methyl-seq/RNbeads_interpretation_visualization.Rmd):
- [dmr-seq.Rmd](https://github.com/stephen-siecinski/social_behavior_in_mice/blob/main/methyl-seq/dmr-seq.Rmd):

## Summary results

<img src="https://github.com/stephen-siecinski/social_behavior_in_mice/blob/main/methyl-seq/dmr.png" width="800">
