# kcyteReg

Analysis scripts assciated with manuscript:

Finnegan A. _et al._"Single-cell transcriptomics reveals spatial and temporal turnover of keratinocyte differentiation regulators."

## Repository Organization
within the ./kcyteReg directory:

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
     + `./codes/corrFuncs.py` and `./codes/run_calcCorr.py`
     +  `./codes/run_clusterCoRegTFcorr.py` 
     + `./codes/run_getClustMapBlocks_inconsistStat.py`
 + Transcriptomic comparison of single-cell population with bulk BCC/SCC samples
     + see  ./BCC_SCC/BCC_SCC_Analysis.ipynb 
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

+ `BCC_SCC` — Transcriptomic comparisong of single-cell stages to DE patterns of BCC, SCC bulk data.
 
__Misc__ 

+ `./raw` — contains files that are starting points for analysis. 

+ `./setsGenes` — Gene sets generated during analysis

## Dependencies
+ Requires python >=3.3, R >= 3.4, git, graphviz (http://www.graphviz.org/) - a pygraphviz dependency which is not automaically installed
+ We recommend installing required python packages into a virtual environment:
	```bash
	## In the project root diretory run:
	python3 -m venv --prompt krUser .  
	source bin/activate       ## use "deactivate" to leave environment
	pip install -r ./install/requirements.txt 
	```
+ To install the appropriate version of MAGIC from github run
	```bash
	 ./install/install_magic0.1.sh 
	```
+ To setup R libraries and evironmental vairables run
	```bash
	./install/setupR.sh local  ## remove "local" if you are not using the virtual environment (not recommended)
	```
+ Set PYTHONPATH by running
	```bash
	source ./.projectRC.sh
	```

+ To run scripts in `./motifAnalysis/` you will need to have FIMO (http://meme-suite.org/doc/fimo.html) and bedtools ( https://bedtools.readthedocs.io/en/latest/content/installation.html) installed.
