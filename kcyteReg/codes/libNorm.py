#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  8 22:35:34 2018

@author: afinneg2
"""
from __future__ import print_function
import numpy as np
import pandas as pd
import argparse

#### Function definitions #####################3
def libNorm_totalExpr(rawCounts, standardSize = 10**4 ,returnLibSizes = False ):
    """
    Normalized express data for each cell to (geneExpr / libSize ) * 10**4
    Dnew_i,j = [D_{i,j} / sum_j(D_{i,j})] * 10**4)
    Inputs:
        rawCounts - dataFrame rows are cells , cols are genes
    Reurns:
        rawCounts  - normalized dataFrame
    """
    libSize = rawCounts.sum(axis=1)
    rawCounts =  rawCounts.divide(libSize , axis = 0) * standardSize
    if returnLibSizes:
        return rawCounts , libSize 
    else:
        return rawCounts
    
#### Parse command-line #################################################
parser = argparse.ArgumentParser(description='Perform library size normalization')
parser.add_argument('--f_i', 
                    help='Expression values. By default cellsn  are rows and columns are genes', required=True)
parser.add_argument("--cellColumns" ,action = "store_true" , 
                    help = "Indicate thta f_i has data for cells in columns instead of rows. \
f_i will be transposed after loading")
parser.add_argument("--standLibSize" , type = float , default = 10.0**4 , 
                    help = "Expression values will be scaled to units per standLibSize")
parser.add_argument('--fo' , 
                    help = "path to which expression values scaled to starndard libSize are written", required= True)
parser.add_argument('--fo_origLibSizes', default = ""  , help  = "path to which original lib sizes are written" )
parser.add_argument('--sep' ,  default = "," , help = "delimiter character for f_i")
parser.add_argument('--debug' , action = "store_true")
args = parser.parse_args() 

#### Load data 
print("Loading data")
data = pd.read_csv(args.f_i , index_col = 0 , sep = args.sep)
print("\tDone")
if args.cellColumns:
    print("transposing loaded maxtrix")
    data = data.transpose()
if args.debug:
    print("the first 10 rows and  columns of loaded data are")
    print(data.iloc[0:10 , 0:10])
#### Perform library size normalizations
print("Normalizing library size to expressoin units per {:.2e}".format(args.standLibSize))
data, libSizes = libNorm_totalExpr(data, standardSize = args.standLibSize, returnLibSizes = True )

#### write outputs
print("Writing results to file")
data.to_csv(args.fo , sep = args.sep, float_format = "%.4f")
if args.fo_origLibSizes:
    libSizes.to_csv(args.fo_origLibSizes, sep = args.sep, float_format = "%.4f"   )
print("\tDone")
