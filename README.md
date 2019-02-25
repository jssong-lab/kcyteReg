# kcyteReg

Analysis scripts assciated with manuscript:

Finnegan A. _et al._"Single-cell transcriptomics reveals spatial and temporal turnover of keratinocyte differentiation regulators."

## Repository Organization

__`codes` — custom codes used for analyis__

 + Single-cell imputation
  + `./codes/run_zinb.R`
  + `./codes/run_MAGIC.py`
 + Dimensionality reduction and  cell clustering
     + `./codes/cluster_spectral.py `and `./codes/run_specCluster.py`
   + `./codes/run_PCA.py`
 + Differential expression
  + `./codes/fantomDE.py`
  + `./codes/run_cluster_DE_subsetCells.R` 
 + Differential motif enrichment
  + `./codes/rankMotifs.py`
  + `./codes/diffEnrich_motifs.py `
 + Identification of gene and TF modules
  + `./codes/corrFuncs.py` and `/codes/run_calcCorr.py`
  +  `./codes/run_clusterCoRegTFcorr.py` 
  + `./codes/run_getClustMapBlocks_inconsistStat.py`
 + Transcriptomic comparison of single-cell population with bulk BCC/SCC samples
  + *\*Fill in\*\* 
 + And other miscellaneous scripts

__Analysis and generation of results__ 

+ `getTFCandidates` — Constructs sets of candidate transcription regulators based on specificity of expression across primary cells

+ `imputeExpr` —scripts for imputation of single cell expression from raw counts

+ `clusterCells` — dimensionality reduction and clustering cells into 8 stages based on imputed expression values 

+ `DE` — One vs rest test for genes differentiatlly expressed in each stage. Tests for genes differentiatlly expressed between the BK and DK state

+ `exprCorr` — Calculation of gene expression correlations across cells in various combinations of stages
+ `motifAnalysis` — Identification of TF binding motifs differentially enriched between SEs unique fo the BK and DK states

+ ` TFexpr_bindingMoitf_assn` — Association between coordinated changes in TF expression and differential binding enrichment between super-enhancers specific to BK and DK states

+ `regNetwork` — Identification of gene and TF modules based on expression correlation 

+ `antiox` — Analysis of expression of genes annotated for antioxidant function across progressive differentiation stages

__Misc__ 

+ `./raw` — contains files that are starting points for analysis. 

+ `./setsGenes` — Gene sets generated during analysis