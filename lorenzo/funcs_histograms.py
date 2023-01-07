# -*- coding: utf-8 -*-
"""
Created on Sun Sep 26 20:07:59 2021

@author: Lorenzo Pica Ciamarra
"""
import os

os.chdir(os.path.dirname(os.path.realpath(__file__)))
os.chdir('..')


import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from itertools import count
from .analysis import get_molecules, get_functionals 
molecule_dictionary = get_molecules()
functional_dictionary = get_functionals()
matplotlib.use('module://ipykernel.pylab.backend_inline')

#%%
def peak_arrays(mol_list, fits):
    """
    

    Parameters
    ----------
    mol_list : TYPE, optional
        List of molecules fitted
    fits : TYPE, optional
        List of fits as output by iterate_mols

    Returns
    -------
    peak_centres_list : TYPE
        Position of the peaks for each molecule.
    peak_ampls_list : TYPE
        Amplitude (area) of the peaks for each molecule.

    """
    peak_centres_list = []
    peak_ampls_list = []
    for i in range(0, len(fits)):
        best = fits[i]
        mol = mol_list[i]
        n_max = best[0]
        model = best[1]
        peak_centres = [model.params['l%d_center' % n] for n in range(1,n_max+1)]
        peak_ampls = [model.params['l%d_amplitude' % n] for n in range(1,n_max+1)]
        zipped = sorted(list(zip(peak_centres, peak_ampls)))
        peak_centres, peak_ampls = list(zip(*zipped))[0][::-1],  list(zip(*zipped))[1][::-1]
        inner_peaks = []
        ampl_inner_peaks = []
        for i in range(0, n_max):
            centre_value = peak_centres[i].value
            ampl_value = peak_ampls[i].value
            inner_peaks.append(centre_value)
            ampl_inner_peaks.append(ampl_value)
        if inner_peaks != []:
             peak_centres_list.append((mol,inner_peaks))
             peak_ampls_list.append((mol, ampl_inner_peaks))
    return peak_centres_list, peak_ampls_list


#Split the region under consideration into several bins. If a molecule has peaks in any given bin, find their average amplitude.
#Sum over the average molecule amplitude into that bin. That's the bin's amplitude.       
def find_bins(low,high, centres, ampls, single_bin = None, k=0):
    w=int((low+high)/300) #width of each bin, determined such that it is related to the wavenumbers (larger wavenumbers, larger bins)
    bounds = list(range(int(low-k),int(high),w)) #Boundaries between the bins
     
    #bins =np.zeros(len(bounds)),[]]
    count_bins  = np.zeros(len(bounds)) 
    ampl_bins = np.zeros(len(bounds)) 
    mol_bins = [[] for x in range(len(count_bins))]
    intervals = []
    for i in range(0, len(centres)): #Take each sublist (index i) in the list of all peaks organised by [[mol0, peaks0], [mol1, peaks1]...] 
        bloc = np.asarray(centres[i][1]) #Take the peaks for mol0
        bloc_ampl = np.asarray(ampls[i][1])
       
        bloc_mol = centres[i][0]#Take the code of mol0
        
        for j in range(0, len(bounds)): #For each molecule, consider a 10-wn (or otherwise calculated) interval (index j)
            if single_bin: #This is only used if we want to check what happens within a certain given bin 
                if type(single_bin) != tuple or len(single_bin) != 2:
                
                    print('If single_bin is given, it must be a tuple of two elements')
                    return None
                elif single_bin[1] <= single_bin[0]:
                    print('The second element of single_bin must be higher than the first')
                    return None
                interval = single_bin
                #print(i)
            elif i == 0: #For each new functional group (run for the first molecule in the list of all molecules with that FG)
                try:
                    interval = (bounds[j], bounds[j+1]) #Consider one of the bins created befor
                except IndexError:
                    interval = (bounds[j], bounds[j]+(bounds[j]-bounds[j-1])) 
                intervals.append(interval)
            else:
                interval = intervals[j]
            in_this_bin = np.where((bloc>interval[0])&(bloc<interval[1]))[0] #Indices of peaks within the bin being considered
            if len(in_this_bin)>0: #If there's any peaks in that bin 
                ampl_this_bin = np.average(np.asarray(bloc_ampl)[in_this_bin]) #Compute the average amplitude of peaks molecule i has in the bin 
                ampl_bins[j] += ampl_this_bin #Add this to the total amplitude of the bin (sum over molecules  i)
                count_bins[j] += 1 #Increase the count of molecules with peaks in that bin by one 
                if bloc_mol not in mol_bins[j]:
                    mol_bins[j].append(bloc_mol) #Record this molecule has a func in that bin 
            if single_bin:
                break #There's only one interval, so there are no other intervals to iterate through
    if single_bin:
        return ampl_bins[0], count_bins[0], mol_bins[0], single_bin 
    zipped = sorted(list(zip(ampl_bins, count_bins, count(), mol_bins, intervals)), key=lambda x: (x[0] is None, x[0]))
    sorted_ampl, sorted_count, sorted_mols, sorted_intervals = list(zip(*zipped))[0][::-1], list(zip(*zipped))[1][::-1], list(zip(*zipped))[3][::-1],  list(zip(*zipped))[4][::-1]  #np.sort(count_bins)[::-1]
    #Return the amplitude of each bin, and the respective counts of how many molecules have peaks there, their names, and the bounds
   # print(np.arrsorted_ampl)
    #print()
    #print(sorted_count)
    return sorted_ampl, sorted_count, sorted_mols, sorted_intervals 





#Check whether there's some degree of degeneracy: do all/some molecules with peaks in a bin also share other functionals 
#with nearby features, which might explain the fact they all share peaks in one bin?
def check_common_funcs(mols, interval, func_tested= '[H]C(=O)[!#1]', verbose=False): #Runs on a single bin
    mol_watchlist = []
    for mol in mols:
        funcs = molecule_dictionary.get(mol)
        mol_fg_codes = []
        for func_tuple in funcs:
            func_code = func_tuple[0]
            func_obj = functional_dictionary[func_code]
            mol_fg_codes.append(func_code)
            if func_code != func_tested:
                #print('check 2')

                for symmetry in func_obj.averageSymmetries():
                    for line in symmetry.properties:
                        fg_min = line.low
                        fg_max = line.high
                        #print(fg_min, fg_max)
                        if (fg_min < interval[1] and fg_max > interval[0]) or (fg_max < interval[1] and fg_max > interval[0]):# Why not considering also only minimum within bounds??
                            if mol not in mol_watchlist:
                                mol_watchlist.append(mol)
    #Assign a "degeneracy score" to the bin given by the proportion of molecules with peaks in that bin that also share other functionals
    #expect to appear in that bin
    try: 
        degeneracy = len(mol_watchlist)/len(mols) 
    except ZeroDivisionError:
        degeneracy=0
    if verbose:
        print('{0} out of {1} molecules in this bin have an additional functional with features here'.format(len(mol_watchlist), len(mols)))
    return degeneracy, mol_watchlist, mols 


#Apply a correction to the bins' amplitudes proportional to their degeneracy score
def multiply_degeneracy(mols, ampls, intervals, func):
    reduced_ampls = np.array([])
    for i in range(0,len(intervals)):
        #print(i)
        degeneracy, mol_watchlist, mols_bin = check_common_funcs(mols = mols[i], interval = intervals[i], func_tested = func)
        reduced_ampl = ampls[i]*(1-degeneracy)
        reduced_ampls = np.append(reduced_ampls, reduced_ampl)  
    return reduced_ampls


def setup_bins(low, high, mols, fits,func = '[H]C(=O)[!#1]', multiply = True, eliminate=False ):
    #print('00 \n',mols)
    centres, ampls = peak_arrays(mol_list=mols,fits=fits)
    sorted_ampl, sorted_count, sorted_mols, sorted_intervals = find_bins(low=low, high=high, centres=centres, ampls=ampls)
   
    ###ADDED TEMP. 13/07/22:##################
    ##########################################
    if type(sorted_count) == tuple:
        sorted_ampl = tuple(np.asarray(sorted_ampl)*np.asarray(sorted_count)) ###
    elif sorted_count*0 == 0: #it means it's a number!
        sorted_ampl = sorted_ampl * sorted_count
    else:
        print('Unable to multiply the original amplitude by the count')
    
    int_avg = np.asarray([(x[0]+x[1])/2 for x in sorted_intervals])
    multiplied_ampls = None
    eliminated_ampls = []
    eliminated_counts = []
    eliminated_mols = [[] for x in range(len(sorted_intervals))]
    if multiply == True:     #Apply the correction proportional to the degeneracy score
        multiplied_ampls = multiply_degeneracy(sorted_mols, sorted_ampl, sorted_intervals, func)
    if eliminate == True: #Simply disregard all molecules with degeneracies when calculating the cumulative amplitudes
        for i in range(0, len(sorted_intervals)):
            degeneracy, mol_watchlist, mols_bin = check_common_funcs(mols = sorted_mols[i], interval = sorted_intervals[i], func_tested = func)
            set_not_deg_mols = set(mols_bin)-set(mol_watchlist)
            ind_not_deg_mols = [mols.index(x) for x in set_not_deg_mols]
            not_deg_mols = [mols[x] for x in ind_not_deg_mols] #Need to do this because the first not_deg_mols is an unordered set
            not_deg_mols_fits = [fits[x] for x in ind_not_deg_mols]
            not_deg_centres, not_deg_ampls = peak_arrays(mol_list=not_deg_mols,fits=not_deg_mols_fits)
          
            elim_ampl, elim_count, elim_mols, interval = find_bins(low=low, high=high, centres=not_deg_centres, ampls=not_deg_ampls, single_bin = sorted_intervals[i])
            ###################################################################
            
            if type(elim_count) == tuple:
                elim_ampl = tuple(np.asarray(elim_ampl)*np.asarray(elim_count)) ###
            elif elim_count*0 == 0: #it means it's a number!
                elim_ampl = elim_ampl * elim_count
            else:
                print('Unable to multiply the eliminated amplitude by the count')
            ###################################################################
            eliminated_ampls.append(elim_ampl)
            eliminated_counts.append(elim_count)
            eliminated_mols[i].append(elim_mols)
    if multiply ==True and eliminate ==True:
        return (sorted_ampl, sorted_count, sorted_mols, sorted_intervals, int_avg), multiplied_ampls, (eliminated_ampls, eliminated_counts, eliminated_mols, sorted_intervals, int_avg)
    elif multiply ==True and eliminate ==False:
        return (sorted_ampl, sorted_count, sorted_mols, sorted_intervals,int_avg), multiplied_ampls
    elif multiply ==False and eliminate ==True:
        return (sorted_ampl, sorted_count, sorted_mols, sorted_intervals,int_avg), (eliminated_ampls, eliminated_counts, eliminated_mols, sorted_intervals, int_avg)
    elif multiply ==False  and eliminate ==False:
        return sorted_ampl, sorted_count, sorted_mols, sorted_intervals, int_avg


#Mostly just plotting and saving the results into files 
def output_histograms(mols, fits, low,high,funcname= '[H]C(=O)[!#1]', plot_hists = False, save_hist_plots = False):
    sorted_ampl, sorted_count, sorted_mols, sorted_intervals, int_avg = setup_bins(low=low, high=high, mols=mols, fits=fits, multiply=False,func=funcname)
    orig, multiplied_ampls = setup_bins(low=low, high=high, mols=mols, fits=fits, multiply=True,func=funcname) 
    orig, elim = setup_bins(low=low, high=high, mols=mols, fits=fits, multiply=False, eliminate=True,func=funcname)
    corrected_funcname = funcname.replace(':','_').replace('\\','£').replace('/','€') #Replace forbidden characters in the filename
    if not os.path.exists(r'..\func_group_stats\{0}'.format(corrected_funcname)):
        os.makedirs(r'..\func_group_stats\{0}'.format(corrected_funcname))
    filename = r'..\func_group_stats\{0}\{1}-{2}'.format(corrected_funcname, int(low), int(high))
    filecontent = np.column_stack((int_avg, sorted_ampl, multiplied_ampls, elim[0], sorted_count, elim[1]))
    np.savetxt('{0}.txt'.format(filename), filecontent, header = '{0}: wn, original ampl, multiplied ampl, eliminated ampl, all count, eliminated count'.format(funcname))#reduced_ampls = multiply_degeneracy(mols=sorted_mols,ampls = sorted_ampl,  intervals = sorted_intervals, func='[H]N([H])[!#1]' )
    if plot_hists ==True:   
        plt.figure()
        plt.title('Original data for {0}'.format(funcname))
        plt.xlabel('Wavenumbers')
        plt.ylabel('Cumulative amplitude')
        plt.bar(int_avg, sorted_ampl, width = 15)
        if save_hist_plots == True:
            plt.savefig('{0}_original.pdf'.format(filename), bbox_inches='tight')
        plt.figure()
        plt.bar(int_avg, multiplied_ampls, width=15)
        plt.title('Data accounting for degeneracy degree for {0}'.format(funcname))
        plt.xlabel('Wavenumbers')
        plt.ylabel('Cumulative amplitude')
        if save_hist_plots == True:
            plt.savefig('{0}_multiplied.pdf'.format(filename), bbox_inches='tight')
        
        plt.figure()
        plt.title('Data deleting degenerate molecules for {0}'.format(funcname))
        plt.xlabel('Wavenumbers')
        plt.ylabel('Cumulative amplitude')
        plt.bar(int_avg, elim[0], width=15)
        if save_hist_plots == True:
            plt.savefig('{0}_eliminated.pdf'.format(filename), bbox_inches='tight')
    return filecontent 