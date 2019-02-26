#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 14:22:37 2018

@author: afinneg2
"""

from __future__ import  print_function

import numpy as np
from scipy.linalg import eigh
from scipy.sparse import issparse, csr_matrix, find
from collections import  OrderedDict
import time

## sklearn
from sklearn.cluster import  KMeans
from sklearn.neighbors import kneighbors_graph
from sklearn.neighbors import NearestNeighbors

### Helper functions ##################################
def calcLap(W, norm = "None"):
    """
    
    Inputs:
        W - square matrix of similarities
        norm - "None" -> comute unnormalized laplacian
    Returns
         Lap , degrees
    """
    if norm == "None":
        diag = np.diag(W)
        degrees = np.sum(W , axis=1)
        Lap = -W
        np.fill_diagonal(Lap,  degrees - diag )
    else:
        raise(NotImplementedError)
    return Lap , degrees

### Versions of spectral clustering
def spectralClust_ShiMalik(W , nClust, eigvals_max = 10,
                          kmeansKwargs = {'init' : "k-means++", "n_init" : 100 , } ):
    """
    W - square matrix of similarities
    eigvals_max - eigen values 
    """
    
    Lap, degrees =  calcLap(W, norm = "None")
    evals , evecs= eigh(a = Lap , b = np.diag(degrees),  
                         eigvals = (0, eigvals_max) , type = 1 ,
                        overwrite_a = True , overwrite_b = True)
    embedding = evecs[:,0:nClust]
    
    ## run k-means
    kmeans = KMeans(n_clusters = nClust , **kmeansKwargs).fit(embedding)
    return evals ,evecs ,  kmeans.labels_,  kmeans.cluster_centers_ , degrees 

##### Approximate spectral clustering ####################
def kasp(X , nClust, alpha , kmeans_frac= None , 
         kmeans1_kwargs = {"n_init" : 20 , "max_iter" : 300, 'njobs' : -1} ,
        kmeans2_maxIter = 300, simFunc = None, simFunc_kwargs = {}):
    """
    Inputs:
        X - nsamples , nFeatures (ndarray)
        nClust - (int)
        alpha - sample reduction factor for spectral clustering (float) 
        kmeans_frac - None or float. 
                        If None run kmeans on all samples. 
                        If float
    
    """
    
    N , m= X.shape
    n_spectral = np.int(np.floor( N / float(alpha) )) ## number of data points used for spectral cluster
    runTimes = OrderedDict([])
    ### Run kmeans
    if  kmeans_frac is not None:
        ### run 2 rounds of kmeans (as in original kasp)
        n_kmeans1 = np.int(kmeans_frac*N)
        if n_kmeans1 < n_spectral:
            raise(ValueError , "fewer data points than clusters for kmeans1")
        kmeans1_idxs = np.random.choice(N, size = n_kmeans1 , replace = False )
        ## kmeans1
        print("started kmeans1")
        ti = time.time()
        kmeans1 =  KMeans(n_clusters = n_spectral , **kmeans1_kwargs).fit(X[kmeans1_idxs , :])
        tf = time.time()
        runTimes["kmeans1 (s)"] = tf - ti
        ## kmeans2
        print("started kmeans2" )
        ti = time.time()
        kmeans2 = KMeans(n_clusters = n_spectral, 
                         init = kmeans1.cluster_centers_ ,
                          n_init = 1, max_iter = kmeans2_maxIter).fit(X)
        tf = time.time()
        runTimes["kmeans2 (s)"] = tf - ti
        
        X_km =  kmeans2.cluster_centers_
        km_labels = kmeans2.labels_
        
    else:
        ## run 1 round kmeans 
        print("started kmeans1")
        ti = time.time()
        kmeans1 =  KMeans(n_clusters = n_spectral , **kmeans1_kwargs).fit(X)
        tf = time.time()
        runTimes["kmeans1 (s)"] = tf - ti
        
        runTimes["kmeans2 (s)"] = None
        
        X_km=  kmeans1.cluster_centers_
        km_labels = kmeans1.labels_
        
    ## construct similarity matrix for X_kmeans 
    print("started constructing similarity" )
    ti = time.time()
    W = simFunc(X_km, **simFunc_kwargs)
    tf = time.time()
    runTimes["construct similarity (s)"] = tf - ti
    
    ## run spectral cluster 
    print("started sepectral cluster" )
    ti = time.time()
    evals, evecs, spec_labels_Xkm, _ , _ =  spectralClust_ShiMalik( W , nClust, eigvals_max = nClust + 10 )
    tf = time.time()
    runTimes["spectral cluster (s)"] = tf - ti
    ## assign class labels to X based on spec_labels_Xkm (spectral cluster labels of Xkm)
    spec_labels = spec_labels_Xkm[ km_labels]
    ### Returning too much information
    return spec_labels, evals, evecs, spec_labels_Xkm, X_km, km_labels , runTimes

### Functions constructing similarity matrices
def knn_similarity(X, n_neighbors = 10 ):
    
    connectivity = kneighbors_graph(X , n_neighbors=n_neighbors, include_self=False)
    # make connectivity symmetric
    connectivity = 0.5 * (connectivity + connectivity.T)
    connectivity = connectivity.toarray()
    return  connectivity     

def compute_magicAffinity(data, k=10, epsilon=1, distance_metric='euclidean', ka=4,
                         returnSparse = False):

    N = data.shape[0]

    # Nearest neighbors
    print('Computing distances')
    nbrs = NearestNeighbors(n_neighbors=k, metric=distance_metric).fit(data)
    distances, indices = nbrs.kneighbors(data)

    if ka > 0:
        print('Autotuning distances')
        for j in reversed(range(N)):
            temp = sorted(distances[j])
            lMaxTempIdxs = min(ka, len(temp))
            if lMaxTempIdxs == 0 or temp[lMaxTempIdxs] == 0:
                distances[j] = 0
            else:
                distances[j] = np.divide(distances[j], temp[lMaxTempIdxs])

    # Adjacency matrix
    print('Computing kernel')
    rows = np.zeros(N * k, dtype=np.int32)
    cols = np.zeros(N * k, dtype=np.int32)
    dists = np.zeros(N * k)
    location = 0
    for i in range(N):
        inds = range(location, location + k)
        rows[inds] = indices[i, :]
        cols[inds] = i
        dists[inds] = distances[i, :]
        location += k
    if epsilon > 0:
        W = csr_matrix( (dists, (rows, cols)), shape=[N, N] )
    else:
        W = csr_matrix( (np.ones(dists.shape), (rows, cols)), shape=[N, N] )

    # Symmetrize W
    W = W + W.T

    if epsilon > 0:
        # Convert to affinity (with selfloops)
        rows, cols, dists = find(W)
        rows = np.append(rows, range(N))
        cols = np.append(cols, range(N))
        dists = np.append(dists/(epsilon ** 2), np.zeros(N))
        W = csr_matrix( (np.exp(-dists), (rows, cols)), shape=[N, N] )
    
    if not returnSparse:
        W =W.toarray()
    
    return W





