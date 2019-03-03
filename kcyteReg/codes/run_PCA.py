#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  5 16:37:46 2018

@author: afinneg2
"""
from __future__ import print_function
import argparse
import os
import numpy as np
import pandas as pd
import sys
import pickle
from sklearn.decomposition import PCA
import time
from IOutils import readListFromFile


#### Function definitions #####################################
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

def run_pca(fitData , projData = None , n_components = 100,  
            svd_solver = 'randomized', returnPCAobj = False,
            scaleFeatures = False):
    if scaleFeatures:
        print("Scaling features")
        means = fitData.mean(axis = 0)
        stds = fitData.std(axis = 0 )
        fitData = (fitData - means ) / stds
        if projData is not None:
            projData = (projData - means )/ stds
        
    pca = PCA(n_components=n_components, svd_solver=svd_solver)
    pca.fit(fitData)
    fitData_pc =  pca.transform(fitData)         ## fitData componets in PC basis 
    loadings = pca.components_                  ## components of PCs in feature basis
    varExplained_ratio = pca.explained_variance_ratio_
    if scaleFeatures:
        np.savetxt("means_test.txt" , pca.mean_)
        assert(np.all(np.isclose(pca.mean_  , 0.0 ) ))
    if projData is None:
        if  returnPCAobj:
            return fitData_pc, loadings, varExplained_ratio, pca
        else:
            return fitData_pc, loadings, varExplained_ratio 
    else:
        projData_pc = pca.transform(projData)    ## projData  componets in PC basis 
        if  returnPCAobj:
            return fitData_pc, loadings, varExplained_ratio,  projData_pc , pca
        else:
            return fitData_pc, loadings, varExplained_ratio,  projData_pc
    
#### Parse command-line #################################################
parser = argparse.ArgumentParser(description='Run PCA')
parser.add_argument('--fi_expr', help='Data file. By default rows are observations columns are features', 
                    required=True)
parser.add_argument("--sep" , default = ",")
parser.add_argument('--nPC' , type =int ,default = 40)
parser.add_argument("--randomized" , action = 'store_true')
parser.add_argument("--geneRows" , action = "store_true" , help = "Indicate thta fi has data for cells in \
columns instead of rows. fi will be transposed after loading")
parser.add_argument("--standLibSize" , type = float , default= 0.0 , help = "If standLibSize >= 1  fi_expr data will be scaled to \
expression units per standLibsize (e.g expression units per 10,000)" )
parser.add_argument('--logT_pseudo' , type = float ,default = -1.0 , help = "pseudo count for log2 transform. \
If pseudo count < 0 no log transform is done (default)")
parser.add_argument("--scaleFeatures" , action = 'store_true' , 
                    help = "Mean center vector of observations for each feature and divide by standard deviation. \
Otherwise only mean centering is performed." )
parser.add_argument('--fi_cellNames' ,  default = "" , help = "names of cells to include . All other cells discarded before PCA" )
parser.add_argument('--fi_geneNames' , default = "" ,help = "names of genes to include. All other genes discarded before PCA")
parser.add_argument('--fi_project' , default = "" , help = "subset of cells in fi_cellNames. \
 Any cells provided we be excluded from PCA fitting and then will be projected onto fitted nPC components.")
parser.add_argument('--fo', help= "output file name. Outputs will be fo-PCcomponents , fo-PCloadings , fo-varExplained")
parser.add_argument('--debug' , action = "store_true")
args = parser.parse_args() 
###### Load Data ############################################################

print("Loading data")
### Load used for PCA
if os.path.splitext(args.fi_expr)[-1] == ".csv":
    data = pd.read_csv(args.fi_expr , index_col = 0 , sep = args.sep)
elif  os.path.splitext(args.fi_expr)[-1] == ".pkl":
    data = pd.read_pickle(args.fi_expr)
else:
    raise Exception("file format {} of fi_expr not recognized".format(os.path.splitext(args.fi_expr)[-1] ))
print("\tDone")
if args.geneRows:
    print("\ttransposing loaded maxtrix")
    data = data.transpose()
if args.debug:
    print("the first 10 rows and  columns of loaded data are")
    print(data.iloc[0:10 , 0:10])
if args.standLibSize >= 1:
    print("Standardizing fi_expr to expression units per {:.3e}".format(args.standLibSize))
    data = libNorm_totalExpr(data, standardSize = args.standLibSize ,returnLibSizes = False )
   
### restrict to specified cells if supplied
if args.fi_cellNames:
    cellNames_allowed = readListFromFile(args.fi_cellNames)
else:
    cellNames_allowed = data.index.values
### restrict to specified genes if supplied
if args.fi_geneNames:
    geneNames_allowed =  readListFromFile(args.fi_geneNames)
else:
    geneNames_allowed  = data.columns.values

nCell_old , nGene_old = data.shape
data = data.loc[cellNames_allowed , geneNames_allowed ].copy()

### log transform if specified by pseudoCount >=0
if args.logT_pseudo >=0:
    print("Applying log2 tranform to data with pseuodo count {:.2f}".format(args.logT_pseudo))
    data = np.log2(data + args.logT_pseudo  )

#### separate Cells used for fitting PCA from cells projected on PC components
if args.fi_project:
    projectOnly = readListFromFile(args.fi_project)
    projectMask =  data.index.isin(projectOnly )
    data_fit =data.loc[~projectMask, :  ].copy()
    data_project = data.loc[projectMask, :].copy()
    nCell_fit , nGene_fit = data_fit.shape
    nCell_project , _  = data_project.shape
    print("Print fitting PCA with {:d} cells (samples) and {:d} genes (features).\n\
Projecting  {:d} cells onto the PCs".format( nCell_fit , nGene_fit,  nCell_project ) )
else:
    nCell_total , nGene_total = data.shape
    print("{:d} of {:d} cells and {:d} of {:d} genes used for PCA".format(nCell_total , nCell_old ,
                                                                          nGene_total , nGene_old))
##### Run PCA ##########################################################################################
if args.randomized:
    svd_solver = 'randomized'
else:
    svd_solver = 'auto'

if args.fi_project:
    ti = time.time()
    fitData_pc, loadings, varExplained_ratio, projData_pc, pca_obj= run_pca(data_fit.values , projData = data_project.values ,
                                                                    n_components = args.nPC,  svd_solver = svd_solver,
                                                                      scaleFeatures =args.scaleFeatures, 
                                                                       returnPCAobj  = True) 
    tf = time.time()
    print("PCA calculation took: {:.3f} s".format(tf-ti))
    ## construct output files 
    pcComps = pd.DataFrame(index = data.index.values ,
                           columns=['PC' + str(i) for i in range(1, args.nPC+1)],
                           dtype =fitData_pc.dtype )
    pcComps.loc[~projectMask, : ] = fitData_pc
    pcComps.loc[projectMask, : ] =  projData_pc  
else:
    ti = time.time()
    data_pc, loadings, varExplained_ratio, pca_obj =  run_pca(fitData = data.values , projData = None ,
                                                    n_components = args.nPC,  svd_solver = svd_solver ,
                                                     scaleFeatures =args.scaleFeatures ,  returnPCAobj = True)
    tf = time.time()
    print("PCA calculation took: {:.3f} s".format(tf-ti))
    pcComps = pd.DataFrame(  data_pc ,index = data.index.values ,
                           columns=['PC' + str(i) for i in range(1, args.nPC+1)])
    
pcLoadings = pd.DataFrame( loadings , 
                          index =  ['PC' + str(i) for i in range(1, args.nPC+1) ] ,
                          columns  =  data.columns )
varExplainedSeries = pd.Series(varExplained_ratio , 
                           index = ["PC{:d}".format(i) for i in np.arange(1 , args.nPC +1 )] )                             
        
#### write to output files ######################################################################
print("Writing results to file")
basename , ext = os.path.splitext(args.fo)
fo_PCcomponents = basename + "-PCcomps" + ext
fo_PCloadings = basename + "-PCloadings" +  ext
fo_varExplained = basename + "-varExplained" + ext
fo_pcaObj = basename + "-PCAobj.pkl"

pcComps.to_csv(fo_PCcomponents)
pcLoadings.to_csv(fo_PCloadings)
varExplainedSeries.to_csv(fo_varExplained)

f = open(fo_pcaObj , 'wb')
pickle.dump(pca_obj , f)
f.close()
print("\tDone")
