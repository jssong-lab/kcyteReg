#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import division
from __future__ import print_function

import numpy as np
from scipy.stats import entropy
import pandas as pd
from collections import OrderedDict
import re
import argparse


#### FUNCTION DEFINITIONS #####################################
def parser_motifFile_toArr(fname):
    motifDict = OrderedDict([])
    buffer = []
    motifName = ""

    f = open(fname , 'r' )
    line = f.readline()

    while line:
        if line.startswith("MOTIF"):
            if motifName:
                motifDict[motifName] = np.array(buffer)
            motifName = line.strip("\n").split()[-1]
            buffer = []
            line  = f.readline()
            continue
        match = re.match(r'\s*([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)' , line )
        if match is not None:
            buffer.append( [float(match.group(i)) for  i in range(1,5)]   )
        line  = f.readline()
    if buffer:
        motifDict[motifName] =np.array( buffer )

    f.close()
    return  motifDict 

def read_background(fname):
    """
    File is of the form
    <line1>
    A <float> C <float> G <float> T <float>
      
    """
    f = open(fname  , 'r')
    lines = f.readlines()
    f.close()
    lines = [x.strip("\n") for x in lines]
    background_freq = [float(x.strip() ) for x in lines[1].split()[1::2] ]
    return np.array(background_freq)

def rank_by_KLdiv( motif_names , motif_dict , bkground = np.array(4*[0.25] ) ):
    
    assert np.sum(bkground) == 1.0
    KLdivs = [ ]
    for motif_name in motif_names:
        p_vec =  motif_dict[motif_name ] ## motifLen x 4
        KLdiv =  np.array([entropy(p , bkground )  for p in  p_vec ]) 
        KLdivs.append( (motif_name , KLdiv )  )
    minLen = min( [ len(x[1])  for x in KLdivs ] )
    KLdivs = list(map( lambda x: ((x[0]), np.sum( sorted(x[1] , reverse = True)[0:minLen] ) ) ,  KLdivs ))
    
    KLdivs_sorted = sorted(KLdivs , key = lambda x:  x[1] , reverse = True)
    motif_names_sorted = [ x[0] for x in KLdivs_sorted ] 
    ranks = np.array([  motif_names_sorted.index(x) for x in motif_names]  ) + 1
    result = pd.DataFrame( np.stack([ ranks  , np.array([x[1] for x in KLdivs] )]).transpose() , index = motif_names)
    result.columns = ["rank" , "KLDiv_sum"]
    result_sorted= result.sort_values(by = "rank" )
    result_sorted.loc[:, "rank"]= result_sorted["rank"].astype("int")
    
    return result_sorted

if __name__ == "__main__":
    ### GLOBALS ####################################################
    #motifHeader_to_symb_fname ="motifHeader_to_scRNAsymbol.txt"  ## Moved to command line argument 
    
    #### PARSE CMD LINE ##############################################
    parser = argparse.ArgumentParser(description= "Prioritize multiple motifs corresponding to the same TF" )
    parser.add_argument("--motifs"  , help = "Motif file in MEME format")
    parser.add_argument("--background", help = "2 line file:\
                        Background letter frequencies\nA <freqA> C 0.25 G 0.25 T 0.25")
    parser.add_argument( "--DBs" , default = "all", help =  "Comma separated list with elements from {JASPAR TRANSFAC Hocomoco Jolma2013} or all")
    parser.add_argument("--motifHeader_to_symb")
    parser.add_argument( "--fo" , help = "output file name" )
    
    #example_args= "--motifs motifs_TFsfantom-Klein-misc.meme \
    #--background ./background.uniform.txt \
    #--DBs all \
    #--fo motifs_TFsfantom-Klein-misc.rank"  
    args = parser.parse_args()
    
    ### LOAD DATA #######################################################
    motif_dict= parser_motifFile_toArr(args.motifs)
    background = read_background(args.background)
    
    ## MotifHeader to symbol
    motifHeader_to_symb = pd.read_csv(args.motifHeader_to_symb , header = 0 ,index_col = 0, sep = "\t")
    motifHeader_to_symb["symbol_scRNAseq"] = motifHeader_to_symb.loc[: , "symbol_scRNAseq"].map( lambda x: ",".join(sorted(x.split(","))) )
    motifHeader_to_symb =  motifHeader_to_symb.loc[list( motif_dict.keys() ) , :].copy()
    
    ##### Group by scRNA symbol(s) and rank within group by KL div #################################
    motifHeader_to_symb_bySymb = motifHeader_to_symb.reset_index().groupby(by = "symbol_scRNAseq")
    motifHeader_to_symb_ranked = motifHeader_to_symb_bySymb.apply(
                        lambda x :  rank_by_KLdiv( list(x.loc[: , "motifHeader"]) ,  motif_dict= motif_dict, bkground= background ) )
    
    ## Add a column  indicating database origin
    motifHeader_to_symb_ranked.reset_index(level = 0 ,inplace=True )
    motifHeader_to_symb_ranked.loc[: , "origin"] = motifHeader_to_symb.loc[ list(motifHeader_to_symb_ranked.index) , "origin"  ]
    
    ##### Write Result ################################################################
    
    motifHeader_to_symb_ranked.to_csv( args.fo , index_label= "motifHeader", sep = "\t", float_format='%.6f' )

