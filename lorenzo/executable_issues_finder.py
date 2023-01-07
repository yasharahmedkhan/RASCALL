# -*- coding: utf-8 -*-
"""
Created on Mon Sep 27 20:05:46 2021

@author: Lorenzo Pica Ciamarra
"""
import os
os.chdir(os.path.dirname(os.path.realpath(__file__)))
os.chdir('..')
from rascall.functional_analysis import perform_analysis, extract_features, find_clusters
from rascall.NIST_funcs import find_suitable_funcs
from rascall.plot_histogram import plot_all_hists
import numpy as np
#%%
def plotting(*args):
    if args:
        all_func_dictionary = args[0]
    else:
        print('Loading functional groups...')
        all_func_dictionary  = find_suitable_funcs(threshold=0,max_threshold=np.inf) 
    original = input('Would you like to plot the histograms before accounting for degeneracy? y/[n] \n')
    multiplication = input('Would you like to plot the histograms accounting for degeneracy by multiplication? [y]/n \n')
    elimination = input("Would you like to plot the histograms  accounting for degeneracy by elimination? y/[n] \n")
    if original == 'y':
        original = True
    else:
        original = False
        
    if multiplication == 'n':
        multiplication = False
    else:
        multiplication = True   
        
    if elimination == 'y':
        elimination = True
    else:
        elimination = False   
    plot_all_hists(func_dict = all_func_dictionary, init=original, multipl = multiplication, elimin = elimination)
    print('Your plots are in the "/all_plots" subdirectory within the main RASCALL directory')

#Find the functional groups with enough molecules so that an analysis can be significant:
print("Welcome to the RASCALL issues finder")
skip = input('If you wish to just plot histograms from previous calculations without running any new ones, type "s". \nOtherwise press enter \n')
if skip =='s':
    plotting()
    quit()
    
nmin = int(input("Only consider functional groups present in at least how many molecules? \n"))
nmax = input("If you wish to set an upper bound to limit computation time, type it here. \nOtherwise, press enter. \n")
try:
    nmax = int(nmax)
    if nmax < nmin:
        print('The upper bound must be larger than the lower bound!')
        print('Proceeding with no upper bound.')
except:
    nmax = np.inf
overwrite = input('Would you like to overwrite previous results? y/[n] \n ')

if overwrite == 'n' or overwrite == '':
    overwrite = False
elif overwrite == 'y':
    overwrite = True
else:
    print('Invalid option selected. Proceeding without overwriting.')
print("Loading functional groups...")   
#%%
func_dictionary = find_suitable_funcs(threshold=nmin,max_threshold=nmax) 
all_func_dictionary = find_suitable_funcs(threshold=0,max_threshold=np.inf) 

all_histograms = {}
#%%
#Iterate through each functional group
i=0
tot = len(func_dictionary.keys())
go=input('Will iterate through {0} functional groups: \n{1} \n \nPress enter to start - be aware, this is SLOW!'.format(tot, func_dictionary.keys()))


for func_name in list(func_dictionary.keys()):
    i+=1
    print('')
    print('Functional group {0}, {1} of {2}'.format(func_name,i, tot))
    #Get the attributes of the corresponding functional object 
    #func_name = 'C=C=C'
    fg_features = extract_features(func_name, func_dictionary) 
    #Split up the various features of each object into different groups, each including features close together
    full_bounds = find_clusters(fg_features) #Upper and lower boundaries of each group, and the group's name (number)
    func_histograms = perform_analysis(full_bounds, func_name, overwrite)
    all_histograms[func_name] = func_histograms
print('Your output is in the "func_group_stats" subdirectory within the main RASCALL directory')

hist_plot=input('Would you like to plot all the histograms, including those calculated previously? [y]/n \n')
if hist_plot == 'y' or hist_plot =='':
    plotting(all_func_dictionary)
    
