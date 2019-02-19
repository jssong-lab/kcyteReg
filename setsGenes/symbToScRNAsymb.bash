#!/bin/bash
set -euo pipefail

if [[ $# == 0 ]]
	then
	printf "Usage\nsymbToScRNAsymb.bash <f_in> <f_out>\n"
	exit
fi
fi_symb=$1
fo=$2

python ../codes/intersect_withSynonLookup.py --fi_query $fi_symb --fi_db ../raw/genesAll.scRNA.txt --fi_synonTable "../raw/gene_synonyms_Hsapiens.txt" --fo $fo



