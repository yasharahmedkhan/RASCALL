
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  8 16:13:18 2021

@author: Lorenzo Pica Ciamarra
"""
import os

os.chdir(os.path.dirname(os.path.realpath(__file__)))
os.chdir('..')

from lmfit.models import VoigtModel, ConstantModel  # LorentzianModel

import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('module://ipykernel.pylab.backend_inline')
import scipy.signal as signal
import numpy as np
from rascall.analysis import get_molecules, get_functionals 
from rascall.NIST_spectra import nist_spectrum
from rascall.plot_NIST import NIST_Smile_List
functional_dictionary = get_functionals()
molecule_dictionary = get_molecules()
NIST_data = NIST_Smile_List()
NIST_Smiles = NIST_data[0]
matplotlib.use('module://ipykernel.pylab.backend_inline')
#%%
#Find all individual peaks in the data as "rigorously" defined, which have at least a certain prominence to avoid peaks likely to
#just be due to noise in the data
    
def find_peaks(wn,data):#,n_asked):
    peaks_ind = signal.find_peaks(data, prominence = 0.1*data)[0]
    peaks_prom = signal.find_peaks(data, prominence=0.1*data)[1]['prominences']
    peaks_wn = wn[peaks_ind]
    peaks_left_ind = np.around(signal.peak_widths(data, peaks_ind, rel_height=0.6)[2]).astype(int) 
    peaks_right_ind = np.around(signal.peak_widths(data, peaks_ind, rel_height=0.6)[3]).astype(int)
    peaks_left_wn = wn[peaks_left_ind]
    peaks_right_wn = wn[peaks_right_ind]
    for i in range(0, len(peaks_left_wn)):
        if i > 0:
            if peaks_left_wn[i] < peaks_right_wn[i-1]:
                if (peaks_right_wn[i] - peaks_left_wn[i]) < (peaks_right_wn[i-1] - peaks_left_wn[i-1]):
                    peaks_right_wn[i-1]  = peaks_left_wn[i] 
                else:
                    peaks_left_wn[i] = peaks_right_wn[i-1]  
    peaks_new = peaks_prom*np.pi*10 #Rough conversion from the prominence of the peaks to their area (can't remember where this comes from)
    return (peaks_left_wn, peaks_right_wn, peaks_wn, peaks_new)

# Apply a simple Voigt fit to each peak found in the above function. individually 
def fit_individual_peaks(peaks, wn, data, *extra):
    #Unpack the output of the above function
    peaks_left_wn, peaks_right_wn, peaks_wn, peaks_new = peaks 
    #Initialise some parameters
    n_small =0
    newdata = data
    peaks_left_wn_new = np.array([])
    peaks_right_wn_new = np.array([])
    peaks_wn_new = np.array([])
    peaks_new_new = np.array([])
    fits = []
    regions_wn=[]
    regions_abs = []
    i=0
    
    nmain = len(peaks_wn) #Count how many rigorously identified ("main") peaks there are
    
    while i in range(0, nmain): #For each main peak
        model = VoigtModel()+ConstantModel() #Create a fitting function (Voigt profile over constant background)
        params = model.make_params() 
        region_ind = np.where((wn >= peaks_left_wn[i]) & (wn <= peaks_right_wn[i])) #Delimit the area of wavenumbers to be fit
        if len(region_ind[0])<4: #If the peak is so narrow that it has fewer points than there are free parameters in the fitting function, expunge it from the list
            n_small += 1
            nmain -=1
            peaks_left_wn, peaks_right_wn=  np.delete(peaks_left_wn,i), np.delete(peaks_right_wn,i)
            peaks_wn, peaks_new =  np.delete(peaks_wn, i), np.delete(peaks_new,i)
            continue
        region_wn = wn[region_ind]
        region_abs = newdata[region_ind]
        params['amplitude'].set(value=peaks_new[i]) #Guess the amplitude at the rough measurement found before
        params['center'].set(value=peaks_wn[i]) #Guess the centre point
        params['sigma'].set(value=peaks_right_wn[i]-peaks_left_wn[i]) #Guess the standard deviation as the peak's width 
        params['c'].set(min=0) #Only allow positive backgrounds
        result = model.fit(region_abs, params, x=region_wn) #Perform the fitting 
        fits.append(result) 
        #Same thing as before, without the peaks expunged from the list because they were too narrow
        peaks_left_wn_new = np.append(peaks_left_wn_new,peaks_left_wn[i])
        peaks_right_wn_new = np.append(peaks_right_wn_new,peaks_right_wn[i])
        peaks_wn_new = np.append(peaks_wn_new,peaks_wn[i])
        peaks_new_new = np.append(peaks_new_new,peaks_new[i])
        best = result.eval(x=wn) 
        newdata = newdata-best #Residuals 
        regions_wn.append(region_wn)
        regions_abs.append(region_abs)
        peaks_upd = (peaks_left_wn_new, peaks_right_wn_new,peaks_wn_new, peaks_new_new)
        i+=1
    if extra:
        return  peaks_upd, nmain, newdata, n_small, fits, regions_wn, regions_abs
    return peaks_upd, nmain, newdata, n_small


#Now we aim at finding peaks which were previously embedded into the main ones ("subpeaks")
def find_subpeaks(n, nmain, wn, newdata, data):
    min_height = 0#.05*data
    n_sub = n - nmain #How many subpeaks do we want to look for?
    if n_sub<=0:
         return ([],[]), False
    subpeaks_ind = signal.find_peaks(newdata, height=min_height)[0] #Find the peaks in the residuals
    if n_sub > len(subpeaks_ind):
        n_sub = len(subpeaks_ind) #If more subpeaks are requested than are in the residuals, lower the number of subpeaks we look for
    subpeaks = newdata[subpeaks_ind] * np.pi * 10 #Rough guess at area
    if n_sub == len(subpeaks):
        all_sub=True
    else:
        all_sub = False
    #Following block of code: if n subpeaks out of m (n<m) are requested, only output data for the n highest subpeaks found (WEAK!)
    subpeaks_highest = np.sort(subpeaks)[::-1][:n_sub]
    subpeaks_consider_ind = subpeaks_ind[np.nonzero(
    np.in1d(subpeaks, subpeaks_highest))]
    subpeaks_consider = newdata[subpeaks_consider_ind]
    subpeaks_wn = wn[subpeaks_consider_ind]
    n_sub =  len(subpeaks_consider) 
    return  (subpeaks_wn, subpeaks_consider), all_sub 



def sq_res(model, data):
    sq_res = (model-data)**2
    sum_sq_res = np.sum(sq_res)
    return sum_sq_res
#Given a list of peaks and subpeaks, with guesses as to their centres, widths, heights etc, find the best fit to them 
def fit_peaks_subpeaks(peaks, subpeaks, wn, data, wn_short, data_short, lowlim, highlim, sidepeaks=True):
    peaks_left_wn, peaks_right_wn, peaks_wn, peaks_new = peaks
    nmain = len(peaks_wn)
    subpeaks_wn, subpeaks_consider = subpeaks #No width guesses for subpeaks as they'd be very unreliable
    n_sub = len(subpeaks_wn)
    n = nmain + n_sub #Total number of (main + secondary) peaks 
    model = VoigtModel(prefix='l1_', nan_policy='omit')
    if sidepeaks: #Allow there to be peaks centred outside the considered region
        for i in range(1, n+2): #Create a fitting function made of a superposition of n Voigt profiles + 2 side ones
            model = model + VoigtModel(prefix='l%d_' % (i+1), nan_policy='omit')
    else:
         for i in range(1, n):
            model = model + VoigtModel(prefix='l%d_' % (i+1), nan_policy='omit')
    model = model + ConstantModel(nan_policy='omit') #And add a constant background >=0
    params = model.make_params()
    params['c'].set(min=-0.00001)#,min=0.00, vary=True) #For some reason 0 isn't accepted 
    #Input guesses 
    for i in range(0, nmain):
        params['l%d_amplitude' % (i+1)].set(value = peaks_new[i], min=0)
        params['l%d_center' % (i+1)].set(value = peaks_wn[i], min = peaks_left_wn[i], max = peaks_right_wn[i])
        params['l%d_sigma' % (i+1)].set(value = (peaks_right_wn[i]-peaks_left_wn[i])/2,max=50)
    #Create maximum guesses for the width and amplitudes of subpeaks (no wider/larger than the widest/largest main peak)
    reasonable_sigma_max = np.amax(peaks_right_wn-peaks_left_wn)
    reasonable_ampl_max = np.amax(peaks_new)
    for i in range(nmain, n):
        params['l%d_amplitude' % (i+1)].set(value = subpeaks_consider[i-nmain]*10*np.pi, min=0, max = reasonable_ampl_max)
        params['l%d_center' % (i+1)].set(value = subpeaks_wn[i-nmain])#, max=3000)#, vary=False, min = 0.99*subpeaks_wn[i-nmain], max = 1.01* subpeaks_wn[i-nmain])  # , vary=False)#min=2700, max=2900)
        params['l%d_sigma' % (i+1)].set(max=reasonable_sigma_max)
    if sidepeaks:    
        params['l%d_amplitude' % (n+1)].set(min=0)
        params['l%d_center' % (n+1)].set(max=lowlim)#, vary=False, min = 0.99*subpeaks_wn[i-nmain], max = 1.01* subpeaks_wn[i-nmain])  # , vary=False)#min=2700, max=2900)
        params['l%d_sigma' % (n+1)].set()
        params['l%d_amplitude' % (n+2)].set(min=0)
        params['l%d_center' % (n+2)].set(min=highlim)#, vary=False, min = 0.99*subpeaks_wn[i-nmain], max = 1.01* subpeaks_wn[i-nmain])  # , vary=False)#min=2700, max=2900)
        params['l%d_sigma' % (n+2)].set()            
    result = model.fit(data, params, x=wn) #Perform the fitting 
    res_model = result.eval(x=wn_short)
    residual = sq_res(res_model,  data_short) #Find the sum of square residuals
    return result, residual

#Put together all of the above (finding the peaks, fitting them individually, finding the subpeaks, performing the complete fit) 
#for a desired total number of peaks
def fit_n_peaks(all_wn,all_data,n,lowlim, highlim,plot=None,  sidepeaks=True):
    func_ind = np.where((all_wn >= lowlim) & (all_wn <= highlim))
    func_ind_large = np.where((all_wn >= 0.98*lowlim) & (all_wn <= 1.02*highlim))
    wn = all_wn[func_ind]
    wn_expand = all_wn[func_ind_large]
    
    data = all_data[func_ind]
    padded_data = pad_with_nan(wn_expand,wn,data)
    peaks = find_peaks(wn,data)#,n)
    peaks, nmain, newdata, n_small = fit_individual_peaks(peaks, wn,data)
    subpeaks, all_sub= find_subpeaks(n, nmain, wn, newdata,data)
    result, residual = fit_peaks_subpeaks( peaks, subpeaks, wn_expand, padded_data,wn,data,lowlim,highlim, sidepeaks=sidepeaks)
    n_tot = nmain + len(subpeaks[0])
    print('Fit for {0} peaks completed'.format(n_tot))
    if plot:
        smooth_x = np.linspace(0, 5000, 5000)
        best_fit = result.eval(x=smooth_x)
        ind_comp = result.eval_components(x=smooth_x)
        plt.figure()
        plt.plot(wn, len(wn)*[result.params['c']], ':', label = 'Constant background')
        plt.plot(wn, data, 'k.', label = 'data')
        plt.plot(smooth_x, best_fit, 'b', label = 'best fit')
        for i in range(1, n_tot+1):
            plt.plot(smooth_x, ind_comp['l%d_' % i])
        plt.xlim(wn[0], wn[-1])
        plt.legend()
    return result, residual,wn,data, all_sub, n_tot#, nmain+n_all_sub 

#Run the n-peaks fit for several values of n between a minimum and a maximum
def fit_peaks(all_wn, all_data, lowlim, highlim, nmin=1, nmax=6,  plot_all = False, sidepeaks=True):
    n = nmin
    fits = []
    residuals = []
    peak_n = []
    frac_red=[]
    while n <= nmax:
        fit, residual,wn,data, all_sub, n_tot= fit_n_peaks(all_wn, all_data, n,plot=plot_all, lowlim=lowlim, highlim=highlim, sidepeaks=sidepeaks)
        fits.append(fit)
        if n==nmin:
            frac_red.append(0.1) #Impose that in order for a 2nd fit to be better than the
        #initial one, it has to outperform it by at least 10%
        else:
            frac_red.append(1-(residual)/residuals[-1])
            if frac_red[-1]<frac_red[-2]:
                break
        if n_tot>n:
            n=n_tot    
        residuals.append(residual)
        peak_n.append(int(n))    
        n += 1
        if all_sub ==True:
            print('NO MORE PEAKS OR SUBPEAKS DETECTED')
            break
    frac_reduction = 1-np.insert(np.asarray(residuals[1:])/np.asarray(residuals[:-1]),0,1) #A measure of by how much each fit is better than the previous one
   #Let i be the progressive numbering of the fits, each with n_i peaks, n_i > n_(i-1), and with r_i sum of square residuals.
   #then the frac_reduction = 1 - r_i/r_(i-1) 
    return peak_n, fits, residuals, frac_reduction, wn, data

#!!! WEAKEST PART IN THIS !!!
#Select the best fit. 
#The best fit j will be the first one for which the frac_reduction f_j is such that f_(j+1) <= f_j
def select_best_fit(all_wn,all_data,nmin=1,nmax=6, lowlim=2687.5, highlim=2906.25, mol='O=CO', plot_all=False, plot_final = False, sidepeaks=True):
    peak_n, fits, residuals, frac_reduction,wn, data = fit_peaks(all_wn,all_data,lowlim, highlim,nmin,nmax,  plot_all=plot_all, sidepeaks=sidepeaks)
    best_n, best_fit, best_residuals = peak_n[0], fits[0], residuals[0]
    for i in range(1,len(frac_reduction)):
        if frac_reduction[i] >= frac_reduction[i-1]:
            best_n = peak_n[i]
            best_fit = fits[i]
            best_residuals = residuals[i]
        else:
            break
    print('\n \n The best fit is achieved by {0} peaks \n \n'.format(best_n))
    smooth_x = np.linspace(0, 5000, 5000)
    best = best_fit.eval(x=smooth_x)
    ind_comp = best_fit.eval_components(x=smooth_x)
    if plot_final:
        plt.figure()
        plt.title(mol)
        plt.plot(wn, len(wn)*[best_fit.params['c']], ':', label = 'Constant background')
        plt.plot(wn, data, '.k', label = 'data')
        plt.plot(smooth_x, best, 'b', label = 'best fit')
        for i in range(1, best_n+1):
            plt.plot(smooth_x, ind_comp['l%d_' % i])
        plt.xlim(wn[0], wn[-1])
        plt.ylim(0,1.5*np.amax(data))
        plt.legend()
    return (best_n, best_fit, best_residuals), (peak_n, fits, residuals, frac_reduction)

#Needed to allow the fitting function to consider peaks whose centre is not included in the wavenumbers region being fitted 
def pad_with_nan(wn_expand,wn,data):
    values_missing = sorted(set(wn_expand)-set(wn))
    pad = np.nonzero(np.in1d(wn_expand,values_missing))[0]
    for i in range(0,len(pad)-1):
        if pad[i]-pad[i+1] != -1:
            pad1=pad[0:i+1]
            pad2=pad[i+1:]
    padded_data = np.pad(data, (len(pad1),len(pad2)), constant_values=np.nan)
    return padded_data

#Go through each molecule which has a specified functional group, group of features by group of features, looking for the number of peaks
#that produces the best fit according to the criteria specified in above definitions
def iterate_mols(funcname='[H]C(=O)[!#1]',lowlim=2687.5,highlim=2906.25, sidepeaks=True, nmin=1, nmax=10):
    usable_molecules = []
    best_fits_info = []
    attempted_fits_info = []
    mols = []
    for molecule_code, molecule_functionals in molecule_dictionary.items():
        if any(funcname in s for s in molecule_dictionary.get(molecule_code)):
            if molecule_code in NIST_Smiles:
                    usable_molecules.append(molecule_code)
    i =0
    for molecule_name in usable_molecules:
        i+=1
        all_wn  = nist_spectrum(molecule_name)[0]
        all_data = nist_spectrum(molecule_name)[1]
        print('\n \n FITTING MOLECULE {0} - number {1} of {2}'.format(molecule_name, i, len(usable_molecules)))
        try:
            best, attempted = select_best_fit(all_wn,all_data,nmin=nmin,nmax=nmax, lowlim=lowlim, highlim=highlim, sidepeaks=sidepeaks, mol=molecule_name)
            best_fits_info.append(best)
            attempted_fits_info.append(attempted)
            mols.append(molecule_name)
        except:
            print('FITTING FAILED!')
            continue
    return mols, best_fits_info, attempted_fits_info