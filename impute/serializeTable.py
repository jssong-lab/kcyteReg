#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  5 

@author: afinneg2
"""

import pandas as pd
import argparse
import os
import sys
sys.path.insert(0,"../codes/")
import IOutils

parser = argparse.ArgumentParser(description= "Convert text table so serialzied file with optional filtering of rows/columns" )
parser.add_argument("--fi" , help = ""  )
parser.add_argument("--sep" , default = ",")
parser.add_argument("--no_indexCol", action = "store_true"  )
parser.add_argument("--no_header" , action = "store_true" ) 
parser.add_argument("--rowNames")
parser.add_argument("--colNames")
parser.add_argument("--fo")
args = parser.parse_args()  

###########################################
### LOAD CSV
if (not args.no_indexCol) and (not args.no_header ):
	df = pd.read_csv( args.fi , sep = args.sep , index_col = 0 , header = 0  )
elif (args.no_indexCol) and (not args.no_header):
	df = pd.read_csv( args.fi , sep = args.sep ,header = 0  )
elif ( not args.no_indexCol) and (args.no_header):
	df = pd.read_csv( args.fi , sep = args.sep, index_col = 0 )
else:
	df = pd.read_csv( args.fi , sep = args.sep )

###########################################
## FILTER INDEX AND COLUMNS
if args.rowNames is not None:
    print("Filtering dataFrame index using allowed rowNames file {}".format(args.rowNames))
    index_use = IOutils.readListFromFile( args.rowNames)
    index_intersect = sorted( set( index_use).intersection( set(df.index) ), key = lambda x: index_use.index(x) )
    print("{:d} of {:d} allowed row names are in dataFrame index".format(
                                                        len(index_intersect),
                                                        len(index_use)))
    df = df.loc[index_intersect, :].copy()
if args. colNames is not None:
    print("Filtering dataFrame columns using allowed colNames file {}".format(args.colNames))
    cols_use = IOutils.readListFromFile( args.colNames)
    cols_intersect = sorted( set( cols_use).intersection( set(df.columns) ), key = lambda x: cols_use.index(x) )
    print("{:d} of {:d} allowed column names are in dataFrame columns".format(
                                                        len(cols_intersect),
                                                        len(cols_use)))
    df = df.loc[: ,cols_intersect].copy()
   
###########################################
## WRITE OUTPUT
if args.fo is  None:
	df.to_pickle( os.path.splitext(args.fi)[0] + ".pkl"  )
elif os.path.splitext(args.fo)[-1] == ".pkl":
	df.to_pickle( args.fo )	
elif os.path.splitext(args.fo)[-1] == ".feather": 
	df.reset_index(inplace = True)
	df.to_feather( args.fo  )
else:
	raise Exception("serialization extension not recognized")






