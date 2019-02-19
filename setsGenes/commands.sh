#!/bin/bash 

#####################################################################################
## Write kcyte_TFs candid_TFs kcyte_genes #################################

#########################################################
## Write FANTOM_TFs
## Selected kcte specific
zcat ../getTFCandidates/fantomDE-TF_keratEpiderm-logFCMean-minExpr0-FDR0.05.stats.tsv.gz | tail -n+2  | \
 awk -F '\t' '{if ( ( $4 + 0.0 ) < 0.05  && ( $7 + 0.0 ) > 0.0  ) { print $1 }  }' > tmp.txt
## convert to scRNAseq gene names
./symbToScRNAsymb.bash tmp.txt TFs_fantomDE.scRNA.txt
rm tmp.txt
## Intersect with genes expressed in scRNAseq data
comm -12 <( cat TFs_fantomDE.scRNA.txt | sort  ) <( cat ../raw/genesGeq1pct_allTissue.kcyte.txt | sort ) > TFs_fantom.txt.tmp
rm TFs_fantomDE.scRNA.txt
## Manual curation:
##  + Remove: 
##	+ CDIP1 -> Despite annotation, lit search does no indicate TF function 
##  + Add HES1 HES4
grep -v 'CDIP1'  TFs_fantom.txt.tmp | cat - <( printf "HES1\nHES4" ) | sort >  TFs_fantom.txt
rm TFs_fantom.txt.tmp

#########################################################
## Write FANTOM_Genes
zcat ../getTFCandidates/fantomDE-allGenes_keratEpiderm-logFCMean-minExpr0-FDR0.05.stats.tsv.gz | tail -n+2  | \
 awk -F '\t' '{if ( ( $4 + 0.0 ) < 0.05  && ( $7 + 0.0 ) > 0.0  ) { print $1 }  }' > tmp.txt
## convert to scRNAseq gene names
./symbToScRNAsymb.bash tmp.txt genes_fantomDE.scRNA.txt
rm tmp.txt
## Intersect with genes expressed in scRNAseq data
comm -12 <( cat genes_fantomDE.scRNA.txt | sort  ) <( cat ../raw/genesGeq1pct_allTissue.kcyte.txt | sort ) > genes_fantom.txt
rm genes_fantomDE.scRNA.txt

########################################################
## Write Klein TFs
comm -12 <( cat ../raw/Klein2016_tableS2.txt | sort  ) <( cat ../raw/genesGeq1pct_allTissue.kcyte.txt | sort ) > TFs_Klein.txt

########################################################################
## Write Gene sets for downstream analyis
## Write Kcyte_TFs
cat TFs_Klein.txt  TFs_fantom.txt | sort | uniq > kcyte_TFs.txt
## Write candid_TFs
comm -23 <( cat TFs_fantom.txt | sort ) <( cat  TFs_Klein.txt | sort ) | sort | uniq  >  candid_TFs.txt 
## Write Kcyte_genes.txt
cat Kcyte_TFs.txt genes_fantom.txt | sort | uniq > kcyte_genes.txt











