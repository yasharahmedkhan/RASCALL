# -*- coding: utf-8 -*-
"""
Created on Sun Sep 26 21:02:44 2021

@author: Lorenzo PicaCiamarra
"""
import os

os.chdir(os.path.dirname(os.path.realpath(__file__)))
os.chdir('..')
from rascall.autofit import iterate_mols
from rascall.funcs_histograms import output_histograms
from itertools import count
import numpy as np
#%%

def extract_features(funcname, functional_dict):
    '''
    

    Parameters
    ----------
    funcname : string
        The name of the functional features are to be extracted from.
    functional_dict : dictionary
        The dictionary containing all functionals and their features

    Returns
    -------
    fg_features : list
        #"line" attribute of the RASCALL functional object.

    '''
    functional = functional_dict[funcname][0]
    fg_avg_wn = []
    fg_features = []
    for sym in functional.averageSymmetries():
        for line in sym.properties:
            avg_wn = line.frequency_average()
            fg_features.append(line)
            fg_avg_wn.append(avg_wn)
    zipped = list(zip(fg_avg_wn, count(), fg_features))#Count is here so the sorting never looks at (unsortable) fg_features
    zipped.sort()
    fg_features = list(zip(*zipped))[2]
    return fg_features


def find_clusters(fg_features):
    clusters = {}
    current_cluster = []
    k=0
    fg_avg_wn = np.asarray([line.frequency_average() for line in fg_features])
    fg_max_wn = np.asarray([line.high for line in fg_features])
    fg_dev = fg_max_wn - fg_avg_wn
    for i in range(0, len(fg_avg_wn)):
        if i == 0:
            dist_wn = 0
            current_cluster.append(fg_features[i])
        else:
            dist_wn = fg_avg_wn[i] - fg_avg_wn[i-1] #Distance between the considered feature and the preceding one
        #A maximum distance for two FG to be in the same group, determined by the width of the FG features
        max_tollerance = (fg_dev[i] + fg_dev[i-1])*2 
        if max_tollerance <=60:
            max_tollerance =60
        if dist_wn > max_tollerance: #If the feature is further from the preceding than a specified distance, include it in a different cluster
            k += 1
            current_cluster = []
        current_cluster.append(fg_features[i])
        current_cluster_name = 'Cluster{0}'.format(k)
        #print('Added feature at {0} wn to {1}'.format(fg_avg_wn[i], current_cluster_name))
        clusters[current_cluster_name] = current_cluster
    full_bounds = []
    #Find the boundaries of each cluster, with some tolerance on either side
    for name,cluster in clusters.items():
        n=1
        lower = cluster[0].low - n*(cluster[0].frequency_average() - cluster[0].low)
        upper = cluster[-1].high + n*(cluster[-1].high - cluster[-1].frequency_average())
        width = upper-lower 
        if width < 100:
            adjust = 100-width/2
            lower -= adjust/2
            upper += adjust/2
        bounds = [lower, upper, name]
        full_bounds.append(bounds)
    return full_bounds #Return a list of the boundaries corresponding to each cluster

def perform_analysis(full_bounds, func_code, overwrite = False):
    histograms = {}
    for i in range(0, len(full_bounds)):
        print('Features cluster {0} of {1} for this functional group.'.format(i+1, len(full_bounds)))
        lowlim = full_bounds[i][0]
        highlim = full_bounds[i][1]
        corrected_funcname = func_code.replace(':','_').replace('\\','£').replace('/','€') #Replace forbidden characters in the filename
        if overwrite == False:
            if os.path.exists(r'..\func_group_stats\{0}\{1}-{2}.txt'.format(corrected_funcname,int(lowlim),int(highlim))):
                print('Pre-existing calculations are present in')
                print(r'..\func_group_stats\{0}\{1}-{2}.txt'.format(corrected_funcname,int(lowlim),int(highlim)))
                print()
                print()
                continue
        mols, bests, atts = iterate_mols(funcname = func_code, lowlim=lowlim, highlim = highlim)
        histo_content = output_histograms(mols, bests, lowlim, highlim,funcname=func_code)
        histograms['{0}-{1}'.format(int(lowlim),int(highlim))] = histo_content
    return histograms

