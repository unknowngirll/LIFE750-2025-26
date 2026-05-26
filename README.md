# LIFE750 Advanced Genomic Analysis - Cycle 2 Report

**Student ID:** 201930340 | **Module:** LIFE750 | **Submission:** May 2026

## Overview
This repository contains the code used to characterize a mutant allele of Gene X and evaluate its impact on gene expression using RNA-seq and CUT&RUN data.

## Scripts
*   **`Part1&3-ubuntoCommands.sh`**: Linux pipeline for variant calling (BWA, FreeBayes, BCFtools) and CUT&RUN integration (bedtools).
*   **`Part2-deseq2_analysis.R`**: R script for differential expression analysis and visualization (DESeq2).

## Software Environment
To recreate the required Ubuntu (WSL) environment, run:

```bash
conda create -n life750 -c bioconda -c conda-forge -y bwa samtools bcftools freebayes vcftools bedtools fastqc
conda activate life750
```

## ALL The codes and Commands

1.  Execute the shell script to run the variant calling and binding site overlap analyses: 
    `bash Part1&3-ubuntoCommands.sh`
2.  Open and run the R script line-by-line in RStudio to perform the differential expression analysis: 
    `Part2-deseq2_analysis.R`

## Key Results
*   **Variant Calling:** Identified 3 variants in mutant Gene X (1 stop-gained at pos 413, 1 synonymous, 1 missense).
*   **Differential Expression:** 69 genes significantly dysregulated (33 up, 36 down; *p*adj < 0.01).
*   **CUT&RUN Integration:** 20 of the 69 differentially expressed genes are directly bound by Gene X (15 up, 5 down).

---
**References**
*   Andrews, S. (2010) 'FastQC'. Available at: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
*   Danecek, P. *et al.* (2021) *GigaScience*, 10(2), giab008.
*   Garrison, E. and Marth, G. (2012) *arXiv preprint arXiv:1207.3907*.
*   Li, H. (2013) *arXiv preprint arXiv:1303.3997*.
*   Love, M.I., Huber, W. and Anders, S. (2014) *Genome Biology*, 15(12), p. 550.
*   Quinlan, A.R. and Hall, I.M. (2010) *Bioinformatics*, 26(6), pp. 841-842.
