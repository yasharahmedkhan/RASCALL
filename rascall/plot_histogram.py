# -*- coding: utf-8 -*-

"""
Created on Thu Jul 14 12:56:43 2022

@author: Lorenzo PicaCiamarra
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os
from rascall.functional_analysis import extract_features
from rascall.analysis import get_functionals, get_molecules
os.chdir(os.path.dirname(os.path.realpath(__file__)))
print(os.getcwd())
def corresponding_rascall(func, func_dict, wn):
    fg_features = extract_features(func, func_dict)
    fg_avg_wn = [fg_features[i].frequency_average() for i in range(0,len(fg_features))]
    low_wn = min(wn)-(wn[1]-wn[0])/2
    high_wn = max(wn)+(wn[-1]-wn[-2])/2
    hist_lines = []
    for line in fg_avg_wn:
        if line<=high_wn and line >=low_wn:
            hist_lines.append(line)
        else:
            continue
    return hist_lines
            
def plot_hist(func, file, func_dict, fileformat = 'new', initial = True, multiplied = True, eliminated = False):
    #print(func)
    initial_dir = os.getcwd()
    if fileformat == 'old':
        wn, whole_ampl, mult_ampl, elim_ampl = np.loadtxt(file, unpack=True)
    elif fileformat == 'new':
        wn, whole_ampl, mult_ampl, elim_ampl, whole_count, elim_count  = np.loadtxt(file, unpack=True)
    else:
        print('Choose a file format, either "old" or "new"')
        return
    corrected_funcname = func.replace(':','_').replace('\\','£').replace('/','€') #Replace forbidden characters in the filename
    hist_lines = corresponding_rascall(func, func_dict, wn)
    os.chdir(os.path.dirname(os.path.realpath(__file__)))
    if not os.path.exists(r'../all_plots'):
        os.makedirs(r'../all_plots')
    if not os.path.exists(r'../all_plots/original'):
        os.makedirs(r'../all_plots/original')
    if not os.path.exists(r'../all_plots/elimination'):
        os.makedirs(r'../all_plots/elimination')
    if not os.path.exists(r'../all_plots/multiplication'):
        os.makedirs(r'../all_plots/multiplication')
    os.chdir(r'../all_plots')
    file = file[:-4] #remove the .txt from the range of wavenumbers

    if initial == True:
        plt.figure()
        plt.title('Original data (degeneracy unaccounted for) for {0} \n Based off {1} molecules'.format(func, func_dict[func][1]))
        plt.xlabel('Wavenumbers')
        plt.ylabel('Cumulative amplitude')
        plt.bar(wn, whole_ampl, width = 15)
        plt.vlines(hist_lines, 0, max(whole_ampl), colors='red', linestyles='dashdot', label = 'RASCALL prediction')
        plt.legend(loc='best')
        plt.savefig('./original/{0}_orig_{1}.jpg'.format(corrected_funcname, file), dpi=300)
        plt.close()
    if multiplied == True:
        plt.figure()
        plt.title('Data (after multiplication) for {0} \nBased off {1} molecules'.format(func, func_dict[func][1]))
        plt.xlabel('Wavenumbers')
        plt.ylabel('Cumulative amplitude')
        plt.bar(wn, mult_ampl, width = 15)
        plt.vlines(hist_lines, 0, max(mult_ampl), colors='red', linestyles='dashdot', label = 'RASCALL prediction')
        plt.savefig('./multiplication/{0}_mult_{1}.jpg'.format(corrected_funcname, file), dpi=300)
        plt.legend(loc='best')
        plt.close()
    if eliminated == True:
        plt.figure()
        plt.title('Data (after elimination) for {0} \nBased off {1} molecules (before elimination)'.format(func, func_dict[func][1]))
        plt.xlabel('Wavenumbers')
        plt.ylabel('Cumulative amplitude')
        plt.bar(wn, elim_ampl, width = 15) 
        plt.vlines(hist_lines, 0, max(elim_ampl), colors='red', linestyles='dashdot', label = 'RASCALL prediction')
        plt.legend(loc='best')
        plt.savefig('./elimination/{0}_elim_{1}.jpg'.format(corrected_funcname, file), dpi=300)
        plt.close()
    #os.chdir(os.path.dirname(os.path.realpath(__file__)))
    #os.chdir(r'../func_group_stats/')
    os.chdir(initial_dir)
    return

def plot_all_hists(func_dict, init=False, multipl=True,elimin=True):   
    matplotlib.use('TkAgg')
    os.chdir(r'../func_group_stats/')
    all_funcs = os.listdir()
   # func_dict = get_functionals()
    # for i in func_dict.keys():
    #     func_dict[i] = (func_dict[i], 0)

    for func in all_funcs:
        os.chdir('./{0}'.format(func))
        #print('In func', os.getcwd())
        all_files = os.listdir()
        #print(all_files)
        for file in all_files:
            if file.endswith('.txt'):
                func_name = func.replace('_',':').replace('£','\\').replace('€','/')
                plot_hist(func_name, file, func_dict, fileformat='new', initial=init, multiplied=multipl, eliminated=elimin)
                #print('Done. next: \n', os.getcwd())
            else:
                continue
        print('Plots for {0} completed'.format(func))
        os.chdir('..')        
    #print(os.getcwd())
    

