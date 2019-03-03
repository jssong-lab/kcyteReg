#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  6 20:31:37 2018

@author: afinneg2
"""

from __future__ import print_function
import numpy as np
import pandas as pd
import argparse
import sys
import os
from IOutils import readListFromFile

####### Function definitions #######################################################
def parseSynonFile(fname ,  officialSymCol , synonCol, synonSep = "|", sep ="\t", header = 0 ):
    """
    officialSymCol -int
    synonCol - int
    """
    symbToSynon_df =pd.read_csv(fname , usecols=np.array([officialSymCol , synonCol]) , skiprows =0 , 
                     sep = sep )
    symbToSynon_df.columns = ["Symbols" , "Synonyms"]

    synonToSymbol = {}
    redundantSynonyms = []
    for _ , row in symbToSynon_df.iterrows():
        synonyms_str  = row.loc["Synonyms"]
        symbol = row.loc["Symbols"]
        if synonyms_str !="-":
            synonyms = list(set(synonyms_str.split(synonSep) ) ) ## make sure unique
            for synonym in synonyms:
                if synonym in synonToSymbol.keys():
                    redundantSynonyms.append(synonym )
                synonToSymbol[synonym ] =  symbol
        synonToSymbol[symbol] = symbol
    redundantSynonyms_uniq = np.unique(redundantSynonyms)
    for synonym  in redundantSynonyms_uniq :
        synonToSymbol.pop(synonym )
    
    return  symbToSynon_df, synonToSymbol, redundantSynonyms_uniq 
    
def searchDB_withSynonLookUp( queryList, dbList, synonLookUp ):
    """
    queryList - list of genes
    dbList - list of official symbols
    synonLookUp - dict of form sunonLookup[synonym] -> official symbol
    """
    inDB = []
    inDB_converted = []
    notInDB = []
    
    for query in queryList:
        if query in dbList:
            inDB.append(query)
        elif query in synonLookUp.keys():
            officialSymb = synonLookUp[query]
            if officialSymb in dbList:
                inDB.append(officialSymb)
                inDB_converted.append(query)
            else:
                alternate_synonyms = [x for x in synonLookUp.keys() if synonLookUp[x] == officialSymb]
                foundAlternateMatch = False
                for x in alternate_synonyms:
                    if x in dbList:
                        inDB.append(x)
                        print("found improper symbol {} in db appending to inDB list".format(x))
                        foundAlternateMatch= True
                        break
                if not foundAlternateMatch:
                    notInDB.append(query)              
        else:
            notInDB.append(query)      
    return inDB ,  inDB_converted  , notInDB


if __name__ == "__main__":            
    #### Parse command-line #######################################################################
    
    parser = argparse.ArgumentParser(description="Intersect query gene list with database (db) gene List, using tabel of synonyms to check for differences in nomenclature")
    parser.add_argument("--fi_query" , help = "Query genes list" )
    parser.add_argument("--fi_db" , help ="database genes list (genes in intersection will be described using nomenclature of this list)")
    parser.add_argument("--fi_synonTable" , help = "table with columns corresponding to official gene symbols and synonyms")
    parser.add_argument("--fo" , help = "")
    args = parser.parse_args() 
    
    ##### Load data #################################################################################
    print("Loading database list")
    DB = readListFromFile(args.fi_db)
    print("Loading query list")
    queryList = readListFromFile(args.fi_query)
    print("Parsing synonmy lookup file")
    _,  synonToSymbol, redundantSynonyms  = parseSynonFile(fname  = args.fi_synonTable,
                                                                   officialSymCol=1 , 
                                                                   synonCol=2, synonSep = "|",
                                                                   sep ="\t", header = 0 )
    print("the lookup file contains {} valid official gene symbol and synonym entries \
and {} invalid entries associated with more than one gene symbol".format( len(synonToSymbol) , len(redundantSynonyms)  ))
    
    ###### Perform intersection ####################################################################
    
    inDB, inDB_converted, notInDB = searchDB_withSynonLookUp( queryList = queryList, 
                                                            dbList = DB , 
                                                            synonLookUp = synonToSymbol)
    
    print("Summary of DB Search:\n\tinDB: {}\n\tinDB and converted with table: {}\n\tnot in DB: {}".format(len(inDB) , 
                                                                                                      len(inDB_converted) , 
                                                                                                         len(notInDB)  ))
    ### write result ###################################3
    f = open(args.fo , 'w')
    for name in inDB:
        f.write(name + "\n")
    f.close()
    print("Done")











