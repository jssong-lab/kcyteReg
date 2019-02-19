#!/bin/bash

#####################################################3
## Workflow to generate ./fantomDE-TF_keratEpiderm-logFCMean-minExpr0-FDR0.05.stats.tsv and ./fantomDE-allGenes_keratEpiderm-logFCMean-minExpr0-FDR0.05.stats.tsv
#1. run code in ./write-hg19.cage_peak_phase1and2combined_tpm_allEntrez-sum.ipynb
#2. run code in  ./callDE_FANTOM_TFs.allGenes.ipynb

#####################################################
gzip ./fantomDE-TF_keratEpiderm-logFCMean-minExpr0-FDR0.05.stats.tsv
gzip ./fantomDE-allGenes_keratEpiderm-logFCMean-minExpr0-FDR0.05.stats.tsv


