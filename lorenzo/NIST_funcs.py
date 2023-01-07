# -*- coding: utf-8 -*-
"""
Created on Sun Sep 26 19:39:07 2021

@author: Lorenzo Pica Ciamarra
"""
import os
os.chdir(os.path.dirname(os.path.realpath(__file__)))
os.chdir('..')
from .analysis import get_molecules, get_functionals
from .plot_NIST import NIST_Smile_List
import numpy as np

def count_NIST(funcname, list=False):
    NIST_data = NIST_Smile_List()
    NIST_Smiles = NIST_data[0]
    molecule_dictionary = get_molecules()
    usable_molecules = []
    for molecule_code, molecule_functionals in molecule_dictionary.items():
            if any(funcname in s for s in molecule_dictionary.get(molecule_code)):
                if molecule_code in NIST_Smiles:
                        usable_molecules.append(molecule_code)
    if list == True:
        return len(usable_molecules), usable_molecules
    else:
        return len(usable_molecules)

def find_suitable_funcs(threshold=10, max_threshold = np.inf):
    functional_dictionary = get_functionals()
    suitable_funcs = {}
    for funcname in functional_dictionary.keys():
        n_mols = count_NIST(funcname)
        if n_mols >= threshold and n_mols <= max_threshold:
            suitable_funcs[funcname] = (functional_dictionary[funcname],n_mols)
    return suitable_funcs
