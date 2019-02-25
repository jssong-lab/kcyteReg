#!/bin/bash

###################################################################################
## Scan SEs with motifs 
## FIMO
SE_files=("../raw/NHEKP-K27Ac-SE.namedInterval.txt" "../raw/NHEKD-K27Ac-SE.namedInterval.txt")
fasta="/home/a-m/afinneg2/igenome/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
motifFile="../raw/motifs_kcyteTFs.meme"
fo_arr=( "FIMO_threshE-4_BK_kcyteTFs" "FIMO_threshE-4_DK_kcyteTFs"  )
thresh=0.0001
background="../raw/background_SE.BK_SE.DK_multiline.txt"

for (( i=0; i<${#SE_files[@]}; i++))
	do
	module load BEDTools/2.26.0-IGB-gcc-4.9.4
	bedtools getfasta -name -fi $fasta -bed ${SE_files[i]} | awk '{ if ( $0 !~/^>/ ) $0=toupper($0); print $0 }' > ${fo_arr[i]}.in.fa
	module purge
	module load MEME/5.0.1-IGB-gcc-4.9.4-Perl-5.24.1-Python-2.7.13
	fimo  --verbosity 3 --bfile $background --max-strand --text --thresh $thresh $motifFile ${fo_arr[i]}.in.fa > ${fo_arr[i]}.txt
	rm ${fo_arr[i]}.in.fa
	module purge  	
done

###################################################################################
## Convert to table of moitf hit counts for each motif in each SE
fimoOur_arr=( "FIMO_threshE-4_BK_kcyteTFs.txt" "FIMO_threshE-4_DK_kcyteTFs.txt"  )
thresh=0.0001 
motifIDLookup="../raw/motifHeader_to_scRNAsymbol.txt"

for fimoHitFile in ${fimoOur_arr[@]}
	do
	o_counts=${fimoHitFile%.txt}.counts.tsv
	python ../codes/fimoOut_to_countTable.py -i $fimoHitFile --thresh $thresh --o_counts $o_counts --motifIDLookup  $motifIDLookup 
done
	
###################################################################################
## Identify SEs unique to BK and DK states
SE_files=("../raw/NHEKP-K27Ac-SE.namedInterval.txt" "../raw/NHEKD-K27Ac-SE.namedInterval.txt")
fo_arr=(BK-SE DK-SE  )

module load BEDTools/2.26.0-IGB-gcc-4.9.4
cat ${SE_files[0]}  | cut -f1-6 > ${fo_arr[0]}.all.bed
cat ${SE_files[1]}  | cut -f1-6 > ${fo_arr[1]}.all.bed
bedtools intersect -v -a ${fo_arr[0]}.all.bed -b ${fo_arr[1]}.all.bed | sort -k1,1 -k2,2n >${fo_arr[0]}.unique.bed 
bedtools intersect -v -a ${fo_arr[1]}.all.bed -b ${fo_arr[0]}.all.bed | sort -k1,1 -k2,2n >${fo_arr[1]}.unique.bed
rm ${fo_arr[0]}.all.bed ${fo_arr[1]}.all.bed

######################################################################################
## Write Identifers of SEs unique to BK or DK state in the same format as fimo output
bed_in="./BK-SE.unique.bed"
fo="./BK-SE.unique.names.txt"
cat $bed_in | awk '{printf "%s::%s:%s-%s\n" , $4 , $1 , $2 , $3 }' | sort | uniq > $fo 

bed_in="./DK-SE.unique.bed"
fo="./DK-SE.unique.names.txt"
cat $bed_in | awk '{printf "%s::%s:%s-%s\n" , $4 , $1 , $2 , $3 }' | sort | uniq > $fo

#####################################################################################
## Rank motifs for the same TF according by decreasing KL divergence from background
motifs="../raw/motifs_kcyteTFs.meme"
background="../raw/background_SE.BK_SE.DK.txt"
motifHeader_to_symb="../raw/motifHeader_to_scRNAsymbol.txt"
fo="./motifs_kcyteTFs.rank.txt"

python ../codes/rankMotifs.py --motifs $motifs --background $background --motifHeader_to_symb $motifHeader_to_symb --fo $fo 

#######################################################################################
## Perform Mann-Whiteny U test for differential enrichment of motifs between unique SE sets

f1="./FIMO_threshE-4_BK_kcyteTFs.counts.lengthNorm.tsv"
f2="./FIMO_threshE-4_DK_kcyteTFs.counts.lengthNorm.tsv"
f1_allowedObs="./BK-SE.unique.names.txt"
f2_allowedObs="./DK-SE.unique.names.txt"
motifRanks="./motifs_kcyteTFs.rank.txt"
o_allMotif="BK-SE.unique_Vs_BK-SE.unique_kcyteTFs.threshE-4.mannWhit.tsv"
o_bestMotif="BK-SE.unique_Vs_BK-SE.unique_kcyteTFs.threshE-4.mannWhit.best.tsv"  
effectMeasure=comonLangEffectSize,differenceOfMedian,zScore,rankBiserialCorrelation 
alpha=0.001

python ../codes/diffEnrich_motifs.py --f1 $f1 --f2 $f2 --f1_allowedObs $f1_allowedObs --f2_allowedObs $f2_allowedObs --motifRanks $motifRanks --effectMeasure $effectMeasure --alpha $alpha --o_allMotif $o_allMotif --o_bestMotif $o_bestMotif --plotEffectMeasure
name=write_$o_allMotif

#######################################################################################
## Write TFs with motifs enriched in BK-SE.unique (negative effect size) and in DK-SE.unique (positive effect size) 
f_i="BK-SE.unique_Vs_BK-SE.unique_kcyteTFs.threshE-4.mannWhit.best.tsv"
f_o="MotifsEnriched_SEs_BK.txt"
tail -n+2 $f_i | awk '{ if ( $7 <  0 ) {print $1  } }' > $f_o

f_o="MotifsEnriched_SEs_DK.txt"
tail -n+2 $f_i | awk '{ if ( $7 >  0 ) {print $1  } }' > $f_o


