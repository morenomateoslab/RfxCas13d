# RfxCas13d
Enhanced RNA-targeting CRISPR-Cas technology in zebrafish

## Introduction

This code was implemented as part of the manuscript [*Enhanced RNA-targeting CRISPR-Cas technology in zebrafish*](https://www.biorxiv.org/content/10.1101/2024.10.08.617220v1). More specifically, it generates the heatmaps from Figures 4B and Extended Data 6, and makes use of the following input files:
```
	- input/guide_info.csv: Efficiency information of guide RNAs used in this study.
	- input/GFP_info.csv: Efficiency information of cell culture data on gfp mRNA from Wessels et al., 2020 (PMID:32518401).
	- input/deepCas13_on_our_gRNAs.csv, input/rnaTargeting_on_our_gRNAs.csv and input/TIGER_on_our_gRNAs.csv: Prediction outputs corresponding to executions of deepCas13, rnaTargeting and TIGER on ourgRNAs.
	- human_GFP_DeepCas13_effic_vs_DeepCas13.csv, human_GFP_tiger_effic_vs_rnaTargeting.csv and human_GFP_tiger_effic_vs_rnaTargeting.csv: Prediction outputs corresponding to executions of deepCas13, rnaTargeting and TIGER on cell culture data on gfp mRNA from Wessels et al., 2020.
```

## Usage

The scripts should be run as follows:

```
Rscript 0_prepare_GFP_correl.R
Rscript 1_corr_heatmaps.R
```
## Requirements
- [R >= 3.6](https://cran.r-project.org/)

### R libraries:
- [pheatmap](https://cran.r-project.org/package=pheatmap)
