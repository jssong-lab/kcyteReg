#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import division
from __future__ import print_function
import numpy as np
import pandas as pd
import os
from collections import OrderedDict

##############################################
## LOAD DATA
def readListFromFile(fi, noEmptyEnd = True):
    """
    Load the lines of file fi as a list
    """
    f = open(fi , 'r' )
    lines_raw = f.readlines()
    f.close()
    lines = [ elem.strip('\n')  for elem in lines_raw ]
    
    if noEmptyEnd:
        while lines[-1]=="":
            lines.pop()
            if len(lines) ==0:
                break
    return lines

def parseStrToDict( stringDict , valueType = "str", pairSep = ";" ):
    keyValue_paris = [tuple(pair.split(":")) for pair in stringDict.split(pairSep ) ]
    if valueType == "int":
        rDict = OrderedDict([ (key, int(value)) for key , value in   keyValue_paris ] )
    elif valueType =="float":
        rDict = OrderedDict([ (key, float(value)) for key , value in   keyValue_paris ] )
    elif valueType == "str":
        rDict = OrderedDict([ (key, value) for key , value in   keyValue_paris ] )
    else:
        raise ValueError("valueType {} not recognized".format(valueType))
    return rDict

def loadDF(fi , read_csv_kws = {"header": 0 , "index_col": 0 , "sep" : ","}):
    
    if os.path.splitext(fi )[-1] in [".csv", ".tsv" , ".txt"]:
        df = pd.read_csv( fi, **read_csv_kws )
    elif os.path.splitext(fi)[-1] == ".pkl":
        df = pd.read_pickle(fi)
    elif os.path.splitext(fi)[-1] == ".feather":
        df = pd.read_feather(fi)
        try:
            index_col = read_csv_kws["index_col"]
            df.set_index( list(df.columns)[index_col] , inplcae = True  )
        except KeyError:
            pass
    return df

def loadCellData( fnameDict = OrderedDict([])  ):
    """
    fnameDict  - values are file names allowed keys are ["PCcomps" , "expr"]. Each file is loaded as
                dataFrame and frames are joined along axis 1 using method "inner".
    """
    dfDict = OrderedDict([])
    print("loading Files")
    for key in fnameDict.keys():
        print("\tloading file {}".format(fnameDict[key]))
        if key == "pcComps":
             dfDict[key] = loadDF(fnameDict[key] )
        elif key.lower() == "expr":
            dfDict[key] = loadDF(fnameDict[key] )
        else:
            raise Exception( "key {} not recognized".format(key))
    if len(dfDict[key].keys()) <= 1:
        df_out = pd.concat(list(dfDict.values()), axis =1 , join = "inner", verify_integrity = True)
    else:
        joinKeysDict = {"expr": "expr" , "pcComps" : "rowData"}
        df_out = pd.concat(list(dfDict.values()), axis =1 , join = "inner",
                           keys = [joinKeysDict[key] for key in dfDict.keys()] , verify_integrity = True)
    ### Give summary of join
    print("Summary of Join:")
    for key in dfDict:
        df = dfDict[key]
        print( "\t{} : {} of {} in join".format( key ,
                                                len( set(df.index).intersection(df_out.index) ) ,
                                                len(df.index) ))
    return  df_out 

