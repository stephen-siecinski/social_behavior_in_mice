# Socially Divergent Phenotypes in Inbred Mouse Models

This project compared transcriptional and epigenetic profiles of the amygdala and hippocampus between C57BL/6J and C58/J inbred mouse strains and between high-sociability and low-sociability C58/J mice. RNA and DNA were isolated simultaneously from tissue punches from the respective brain regions and libraries were generated using the Illumina TrueSeq mRNA LP and Agilent SureSelect XT target-enrichment kits. 

The analyses are dived into three primary directories corresponding to distinct pipelines used to interpret the resulting data
- [Cellular Deconvolution](https://github.com/stephen-siecinski/social_behavior_in_mice/tree/main/cellular_deconvolution)
  - Utilizes a recently developed algorithm (MuSiC) to estimate proportions of annotated cell types within bulk-seq data using reference single cell RNA-seq datasets. 
- [RNA-seq](https://github.com/stephen-siecinski/social_behavior_in_mice/tree/main/rna-seq)
- [Methyl-seq](https://github.com/stephen-siecinski/social_behavior_in_mice/tree/main/methyl-seq)
