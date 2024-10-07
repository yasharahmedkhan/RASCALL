# -*- coding: utf-8 -*-
"""
ATMOS 1 - plotting the spectral features of individual functionals within molecules
"""

import pickle as pickle 
import re

from .functional_parser import Functional_Parser
from .molecule_parser import Molecule_Parser
from .molecule_filter import Molecule_Filter
from .plotter import Plotter
from .catalogue import Catalogue
from . import get_file
from .atmospheric_windows import ATMOSPHERIC_WINDOWS
from .molecule_filter import Molecule_Filter

from . import NIST_spectra

from .plot_NIST import NIST_Smile_List

def get_functionals(filename='functionals_formatted_eye_edit.csv'):
    # Load Functionals
    # Example data
    # COC C-O-C sbend 2500 2720 weak
    # COC C-O-C abend 2800 2920 strong
    with open(get_file('functionals.csv'), 'rU') as fhl:
        return Functional_Parser().functional_dictionary_for(fhl.readlines())

#print 'Total number of unique functionals', len(functional_dictionary)

def get_molecules(filename=get_file('dictfunct.p')):
    # Load Molecules
    #Molecule dictionary sample [('C(C1)(C1F)(CC)', [('[H]C([H])(C)C', 2), ('[H]C([H])([!#1])[!#1]', 2),('[H]C([H])([H])C', 1), ('[H]C([H])([H])[!#1]', 1)]),...]
    return pickle.load(open(filename, "rb"))


#print 'Molecule dictionary size', len(molecule_dictionary.items()), 'with sample:', molecule_dictionary.items()[:8]
#print 'Functionals for molecule C(C)NCC(O)', molecule_dictionary.get('C(C)NCC(O)')

def get_plotables(filename='plotable_molecules'):
    # looks through all of the plottable molecules
    plotables = []
    plotable_molecules = open('plotable_molecules', "r")
    for line in plotable_molecules:
        columns = line.strip().split()
        plotables.append(columns[0])
    return plotables



#Makes a list of all hydrocarbons
def get_hydrocarbons(plotables, molecule_dictionary):
    not_hydrocarbons = 0
    hydrocarbons = 0
    hydrocarbons_list = []
    hydrocarbons_in_NIST = []

    for molecule_code, molecule_functionals in molecule_dictionary.items():
        if any('O' in s for s in molecule_code) or any('P' in s for s in molecule_code) or any('N' in s for s in molecule_code) or any('S' in s for s in molecule_code) or any('F' in s for s in molecule_code) or any('I' in s for s in molecule_code) or any('B' in s for s in molecule_code) or any('l' in s for s in molecule_code):
            not_hydrocarbons = not_hydrocarbons + 1
        else:
            #print (hydrocarbons + 1), molecule_code, ' with ', len(molecule_dictionary.get(molecule_code)), ' functionals'
            hydrocarbons_list.append(molecule_code)
            hydrocarbons = hydrocarbons + 1
            if molecule_code in plotables:
                hydrocarbons_in_NIST.append(molecule_code)
    print ('Hydrocarbons:', len(hydrocarbons_list),'Hydrocarbons in NIST:', len(hydrocarbons_in_NIST))
    return hydrocarbons_in_NIST, not_hydrocarbons

def get_halocarbons(plotables, molecule_dictionary):
    not_hydrocarbons = 0
    halocarbons = 0
    halocarbons_list = []
    halocarbons_in_NIST = []

    for molecule_code, molecule_functionals in molecule_dictionary.items():
        if any('O' in s for s in molecule_code) or any('P' in s for s in molecule_code) or any('N' in s for s in molecule_code) or any('S' in s for s in molecule_code):
            not_hydrocarbons = not_hydrocarbons + 1
        elif any('F' in s for s in molecule_code) or any('I' in s for s in molecule_code) or any('B' in s for s in molecule_code) or any('l' in s for s in molecule_code):
            #print (hydrocarbons + 1), molecule_code, ' with ', len(molecule_dictionary.get(molecule_code)), ' functionals'
            halocarbons_list.append(molecule_code)
            halocarbons = halocarbons + 1
            if molecule_code in plotables:
                halocarbons_in_NIST.append(molecule_code)
    print ('Halocarbons:', len(halocarbons_list),'Halocarbons in NIST:', len(halocarbons_in_NIST))
    return halocarbons_in_NIST, not_hydrocarbons

def filter_all(mol):
    return True

def filter_halo(mol):
    return any([elem in mol for elem in ['Fl', 'I', 'B', 'Cl']]) \
        and not any([elem in mol for elem in ['O', 'P', 'N', 'S']])

LETTERS = re.compile(r'\w')
def filter_hydro(mol):
    return set(LETTERS.findall(mol)) == {'C'}


MOL_FILTERS = {
    "all": filter_all,
    "halo": filter_halo,
    "hydro": filter_hydro,
}



# def plot(functional_to_test=None, molecule_fam="all", single_molecule_to_search_for=None, atmospheric_window=None):
#     results = []
#     plotter = Plotter()
#     results.append(f'Functional to test: {functional_to_test}')

#     NIST_data = NIST_Smile_List()
#     NIST_Smiles = NIST_data[0]
#     functional_dictionary = get_functionals()
#     molecule_dictionary = get_molecules()
#     all_molecule_codes = list(molecule_dictionary.keys())
#     molecules = Molecule_Parser().molecules_for(molecule_dictionary, functional_dictionary)
#     molecule_filter = Molecule_Filter(molecules)
#     results.append(f'Total Number of Functional Groups {len(functional_dictionary)}')

#     if atmospheric_window:
#         if atmospheric_window not in ATMOSPHERIC_WINDOWS:
#             results.append(f"Unknown atmospheric window: {atmospheric_window}")
#             return results
        
#         filtered_molecules = molecule_filter.filter_for_atmospheric_window(atmospheric_window)
#         results.append(f"Number of molecules in {atmospheric_window} atmospheric window: {len(filtered_molecules)}")
#         for molecule_code in filtered_molecules:
#             fig = plotter.create_figure(f"{molecule_code} Spectra", molecule_code)
#             plotter.plot_molecule_band_centers(molecules[molecule_code], fig)
#             if molecule_code in NIST_Smiles:
#                 plotter.plot_NIST_spectrum(molecule_code, fig)
#         return results

#     # Handle `single_molecule_to_search_for` argument:
#     if single_molecule_to_search_for != None:
#         mol = molecules[single_molecule_to_search_for]
#         print(mol)
#         return
    
#     if functional_to_test in ("all", "database"):
#         results.append('Plotting all molecules is not recommended due to the large dataset.')
#         return results    
    

#     # Handle `functional_to_test` argument
#     if functional_to_test == "all" or functional_to_test == "database":
#         molecules_with_test_functional = all_molecule_codes
#         molecules_with_test_functional_in_NIST = NIST_Smiles
#         print ('Number of molecules in RASCALL:', len(molecules_with_test_functional))
#         print ('Number of molecules in NIST:', len(molecules_with_test_functional_in_NIST))

#         if functional_to_test == "database":
#             for molecule_code in molecules_with_test_functional:
#                 counter = counter + 1
#                 if len(molecule_dictionary.get(molecule_code)) > 0:
#                     if counter % 100 == 0 or counter == len(molecule_dictionary): print(counter)
#                     plotter.plot_molecule_band_centers(molecules[molecule_code])
#                 else:
#                     print (molecule_code, 'has no functionals')
#         return
#     elif functional_to_test != None:
#         for molecule_code, molecule_functionals in molecule_dictionary.items():
#             if any(functional_to_test in s for s in molecule_dictionary.get(molecule_code)):
#                 if len(molecule_dictionary.get(molecule_code)) > 0:
#                     molecules_with_test_functional.append(molecule_code)
#                 else:
#                     print (molecule_code, 'has no functionals')
#                 if molecule_code in NIST_Smiles:
#                     molecules_with_test_functional_in_NIST.append(molecule_code)

#         if len(molecules_with_test_functional) == 0:
#             print ('Requested functional group', functional_to_test, 'not in database')
#             return

#         print ('Number of molecules_with_test_functional:', len(molecules_with_test_functional))
#         for molecule in molecules_with_test_functional:
#             print(molecule)
#         print ('Number of molecules_with_test_functional in NIST:', len(molecules_with_test_functional_in_NIST))

#         if len(molecules_with_test_functional_in_NIST) < 400:
#             print ('Molecules with this functional also in NIST:', molecules_with_test_functional_in_NIST)
#         else:
#             print ('Too many molecules with this functional to list')
#         return
   
#     # Handle `molecule_fam` argument
#     molecules_in_family = [
#         code for code in all_molecule_codes if MOL_FILTERS[molecule_fam](code)
#     ]
#     molecules_in_family_in_NIST = [
#         code for code in molecules_in_family if code in NIST_Smiles
#     ]

#     results.append(f"Number of RASCALL molecules in family {molecule_fam}: {len(molecules_in_family)}")
#     results.append(f"Number of NIST molecules in family {molecule_fam}: {len(molecules_in_family_in_NIST)}")

#     for molecule_code in molecules_in_family_in_NIST:
#         fig = plotter.create_figure(f"{molecule_code} Spectra", molecule_code)
#         plotter.plot_molecule_band_centers(molecules[molecule_code], fig)
#         plotter.plot_NIST_spectrum(molecule_code, fig)
    
#     return results

import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def plot(functional_to_test=None, molecule_fam="all", single_molecule_to_search_for=None, atmospheric_window=None):
    results = []
    plotter = Plotter()
    logger.info(f'Functional to test: {functional_to_test}')
    results.append(f'Functional to test: {functional_to_test}')

    NIST_data = NIST_Smile_List()
    NIST_Smiles = NIST_data[0]
    functional_dictionary = get_functionals()
    molecule_dictionary = get_molecules()
    all_molecule_codes = list(molecule_dictionary.keys())
    molecules = Molecule_Parser().molecules_for(molecule_dictionary, functional_dictionary)
    molecule_filter = Molecule_Filter(molecules)
    logger.info(f'Total Number of Functional Groups {len(functional_dictionary)}')
    results.append(f'Total Number of Functional Groups {len(functional_dictionary)}')

    molecules_with_test_functional = []
    molecules_with_test_functional_in_NIST = []

    if atmospheric_window:
        if atmospheric_window not in ATMOSPHERIC_WINDOWS:
            logger.warning(f"Unknown atmospheric window: {atmospheric_window}")
            results.append(f"Unknown atmospheric window: {atmospheric_window}")
            return results
        
        filtered_molecules = molecule_filter.filter_for_atmospheric_window(atmospheric_window)
        logger.info(f"Number of molecules in {atmospheric_window} atmospheric window: {len(filtered_molecules)}")
        results.append(f"Number of molecules in {atmospheric_window} atmospheric window: {len(filtered_molecules)}")
        for molecule_code in filtered_molecules:
            fig = plotter.create_figure(f"{molecule_code} Spectra", molecule_code)
            plotter.plot_molecule_band_centers(molecules[molecule_code], fig)
            if molecule_code in NIST_Smiles:
                plotter.plot_NIST_spectrum(molecule_code, fig)
        return results

    if single_molecule_to_search_for is not None:
        logger.info(f"Searching for single molecule: {single_molecule_to_search_for}")
        results.append(f"Searching for single molecule: {single_molecule_to_search_for}")
        if single_molecule_to_search_for in molecules:
            mol = molecules[single_molecule_to_search_for]
            logger.info(f"Found molecule: {mol}")
            results.append(f"Found molecule: {mol}")
        else:
            logger.warning(f"Molecule not found: {single_molecule_to_search_for}")
            results.append(f"Molecule not found: {single_molecule_to_search_for}")
        return results
    
    if functional_to_test in ("all", "database"):
        logger.warning('Plotting all molecules is not recommended due to the large dataset.')
        results.append('Plotting all molecules is not recommended due to the large dataset.')
        return results    
    
    if functional_to_test is not None:
        for molecule_code, molecule_functionals in molecule_dictionary.items():
            if any(functional_to_test in s for s in molecule_dictionary.get(molecule_code)):
                if len(molecule_dictionary.get(molecule_code)) > 0:
                    molecules_with_test_functional.append(molecule_code)
                else:
                    logger.info(f"{molecule_code} has no functionals")
                if molecule_code in NIST_Smiles:
                    molecules_with_test_functional_in_NIST.append(molecule_code)

        if len(molecules_with_test_functional) == 0:
            logger.warning(f'Requested functional group {functional_to_test} not in database')
            results.append(f'Requested functional group {functional_to_test} not in database')
            return results

        logger.info(f'Number of molecules_with_test_functional: {len(molecules_with_test_functional)}')
        results.append(f'Number of molecules_with_test_functional: {len(molecules_with_test_functional)}')
        logger.info(f'Number of molecules_with_test_functional in NIST: {len(molecules_with_test_functional_in_NIST)}')
        results.append(f'Number of molecules_with_test_functional in NIST: {len(molecules_with_test_functional_in_NIST)}')

        if len(molecules_with_test_functional_in_NIST) < 400:
            logger.info(f'Molecules with this functional also in NIST: {molecules_with_test_functional_in_NIST}')
            results.append(f'Molecules with this functional also in NIST: {molecules_with_test_functional_in_NIST}')
        else:
            logger.info('Too many molecules with this functional to list')
            results.append('Too many molecules with this functional to list')
        return results
   
    # Handle `molecule_fam` argument
    molecules_in_family = [
        code for code in all_molecule_codes if MOL_FILTERS[molecule_fam](code)
    ]
    molecules_in_family_in_NIST = [
        code for code in molecules_in_family if code in NIST_Smiles
    ]

    logger.info(f"Number of RASCALL molecules in family {molecule_fam}: {len(molecules_in_family)}")
    results.append(f"Number of RASCALL molecules in family {molecule_fam}: {len(molecules_in_family)}")
    logger.info(f"Number of NIST molecules in family {molecule_fam}: {len(molecules_in_family_in_NIST)}")
    results.append(f"Number of NIST molecules in family {molecule_fam}: {len(molecules_in_family_in_NIST)}")

    for molecule_code in molecules_in_family_in_NIST:
        fig = plotter.create_figure(f"{molecule_code} Spectra", molecule_code)
        plotter.plot_molecule_band_centers(molecules[molecule_code], fig)
        plotter.plot_NIST_spectrum(molecule_code, fig)
    
    return results




# ATMOSPHERIC_WINDOWS = {
#     "co2": [(420.4, 524.0), (822.4, 922.8), (992.8, 1092.8), (1099.6, 1890.0),
#             (1941.6, 2042.0), (2162.0, 2262.0), (2420.0, 3482.0), (3784.0, 4736.0),
#             (5172.0, 6072.0), (6122.0, 6222.0), (6266.0, 6366.0), (6386.0, 6892.0),
#             (7014.0, 8176.0), (8208.0, 8308.0)],
#     "earth": [(430.4, 530.4), (530.8, 630.8), (705.2, 1335.6), (1344.0, 1444.0),
#               (1585.6, 1685.6), (1816.0, 1916.0), (1928.0, 2284.0), (2404.0, 3504.0),
#               (3514.0, 3614.0), (3976.0, 5208.0), (5490.0, 7116.0), (7132.0, 7232.0)],
#     "methane": [(420.4, 1042.4), (1044.4, 1144.4), (1861.6, 2272.0), (3320.0, 3636.0),
#                 (4754.0, 5082.0), (5158.0, 5258.0), (6268.0, 6610.0), (6612.0, 6712.0),
#                 (7704.0, 8116.0), (8144.0, 8244.0), (9058.0, 11020.0)]
# }

#Sample atmosphere windows, hardcoded. Includes CO2-rich, Earth-like, CH4-rich, and Pre-biotic Earth-like atmospheres 
# 'print co2_windows[0][0]'' gives the first item of the first tuple (low freq of the first window) of a CO2 atmosphere

# #co2 atmosphere
# co2_atmosphere = [(420.4,524.0),(822.4,922.8),(992.8,1092.8),(1099.6,1890.0),(1941.6,2042.0),(2162.0,2262.0),
#               (2420.0,3482.0),(3784.0,4736.0),(5172.0,6072.0),(6122.0,6222.0),(6266.0,6366.0),
#               (6386.0,6892.0),(7014.0,8176.0),(8208.0,8308.0)]

# #Earth atmosphere
# earth_atmosphere = [(430.4,530.4),(530.8,630.8),(705.2,1335.6),(1344.0,1444.0),(1585.6,1685.6),
#                     (1816.0,1916.0),(1928.0,2284.0),(2404.0,3504.0),(3514.0,3614.0),
#                     (3976.0,5208.0),(5490.0,7116.0),(7132.0,7232.0)]
# earth2_atmosphere = [(802.8,972.4),(2402.0,2796.0),(4434.0,4806.0),(5626.0,6702.0)] #plus out of range windows between 7540.0-14765.0



# #Methane atmosphere
# methane_atmosphere= [(420.4,1042.4),(1044.4,1144.4),(1861.6,2272.0),(3320.0,3636.0),
#                     (4754.0,5082.0),(5158.0,5258.0),(6268.0,6610.0),(6612.0,6712.0),
#                     (7704.0,8116.0),(8144.0,8244.0),(9058.0,11020.0)]

# methane2 = [(420.4,992.4),(1902.8,2282.0),(3362.0,3534.0),(4894.0,5002.0),(7920.0,8022.0),(9186.0,10900.0),(11545.0,14350.0)]

# #Prebiotic Earth-like windows
# prebio_earth = [(5178.0,5278.0),(6550.0,6866.0),(7630.0,8132.0),(8946.0,11175.0),(11380.0,12315.0)]

# prebio_earth2 = [(420.4, 525.6), (1632.0, 1885.6), (2414.0, 2546.0), (3206.0, 3360.0), (3366.0, 3482.0),
#                  (3934.0, 4034.0), (4630.0, 4736.0), (5166.0, 5316.0), (5320.0, 5558.0), (6374.0, 6892.0),
#                  (7462.0, 8472.0), (8730.0, 15830.0)]

# #features of HCN, approximately. See hcn.agr for spectra at three resolutions
# hcn_regions= [(0.0,100.0),(600.0,820.0),(1340.0,1500.0),(3200.0,3400.0)]
# hcn_strong = [(0.0,100.0),(600.0,820.0)]
# hcn_strong_infrared = [(600.0,820.0)]

# # L, M, N, and Q atmospheric windows
# l_band = [(2710.0, 3115.3)]
# m_band = [(2008.0,2207.5)]
# n_band = [(851.1, 1081.1)]
# lmnq_windows = [(851.1, 1081.1),(2008.0,2207.5),(2710.0, 3115.3)]

# #test_window
# test_window = [(600,950),(1400,1500)]

# # Finds all the strong features whose average frequency is within a window
# atmosphere = methane2
# intensity = 3

# #molecule_filter = Molecule_Filter(molecules)
# #filtered_molecules = molecule_filter.filter_for_region(atmosphere)


# #strength_filtered_molecules = molecule_filter.filter_for_region_and_intensity(atmosphere, intensity)
# #print "our return", len(strength_filtered_molecules)


# #strong_window_molecules = filtered_molecules[1]

                    
# #                for point in window_filtered_list:
# #  #      print point
# #                window_molecules.append(molecule.code)
# #                if point[1] >= 3:
                    
#  #           window_filtered_list = list(filter(lambda x: ((window_low < low_frequency < window_high) and (window_low < high_frequency < window_high)), molecule.high_and_low_frequencies()))
#             #print molecule.high_and_low_frequencies()[0][0],molecule.high_and_low_frequencies()[0][1], molecule.average_points()[0][0]
# #    low_filtered_list = list(filter(lambda x: 600 < x[0] < 800, molecule.high_and_low_frequencies()[0]))
 

# #and window_low < x[1] < window_high


# #            print point
               
# #    if len(window_filtered_list) > 0:
# #        window_molecules.append(molecule.code)
# #        print points
# #        strong_filtered_list = list(filter(lambda y: y[1] == 3, points))
# #        if len(strong_filtered_list) > 0:
# #            strong_window_molecules.append(molecule_code)
        
# #        print 'Molecule exits in window', molecule.code, (filtered_list)

# #print 'List of Molecules', strong_window_molecules


