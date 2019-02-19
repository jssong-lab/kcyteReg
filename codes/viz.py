#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 17:43:39 2019

@author: afinneg2
"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText

def scatterPlot_colors( xvals , yvals , colors ,ax = None ,
                       xlabel = "", ylabel ="" ,
                       xLower = None, xUpper = None ,
                       yLower = None, yUpper = None ,
                       plotKwargs= {} , figsize = (7.5,3),
                      grid = True,  text = "" , textLoc= 1,
                      title = "" , marker ='o',
                      outlierShiftFactor = 0.01):

    if ax is None:
        fig = plt.figure(figsize = figsize)
        ax = fig.add_subplot(1,1,1)
    markers = np.array([marker]*len(xvals))
    ## handle points outside ax limits
    if (xLower is not None) or (xUpper is not None):
        xvals = xvals.copy()
        xOutlierShift =  (min([ elem for elem in [ xUpper, np.max(xvals)]if elem is not None ]) \
                        - max([ elem for elem in [ xLower, np.min(xvals)] if elem is not None ])) *  outlierShiftFactor
        if xLower is not None:
            lowerMask = xvals < xLower
            xvals = xvals.copy()
            xvals[lowerMask] = xLower  + xOutlierShift
            markers[lowerMask] = "<"
        if xUpper is not None:
            upperMask = xvals > xUpper
            xvals[upperMask] = xUpper -  xOutlierShift
            markers[upperMask] = ">"
    if (yLower is not None) or (yUpper is not None):
        yvals = yvals.copy()
        yOutlierShift = ( min([ elem for elem in [ yUpper, np.max(yvals)]if elem is not None ]) \
                        - max([ elem for elem in [ yLower, np.min(yvals)] if elem is not None ])  )*  outlierShiftFactor
        if yLower is not None:
            lowerMask = yvals < yLower
            yvals[lowerMask] = yLower +  yOutlierShift
            markers[lowerMask] = "v"
        if yUpper is not None:
            upperMask = yvals > yUpper
            yvals[upperMask] =  yUpper -  yOutlierShift
            markers[upperMask] = "^"
    uniqueMarkers = set(markers)
    for um in uniqueMarkers:
        mask = markers ==um
        ax.scatter(xvals[mask] , yvals[mask] , c = colors[mask] , marker = um,  **plotKwargs)

    if xlabel: ax.set_xlabel(xlabel)
    if ylabel: ax.set_ylabel(ylabel)
    if (xLower is not None) or (xUpper is not None): ax.set_xlim(left = xLower , right= xUpper)
    if (yLower is not None) or (yUpper is not None): ax.set_ylim(bottom = yLower ,top = yUpper )
    if title: ax.set_title(title)

    if text:
        anchored_text =AnchoredText(text , loc =textLoc)
        ax.add_artist(anchored_text)

    ax.grid(grid)
    return ax
