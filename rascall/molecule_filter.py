import numpy as np
from .atmospheric_windows import ATMOSPHERIC_WINDOWS

class Molecule_Filter:
    def __init__(self, molecules):
        self.molecules = molecules

    # Existing methods
    def filter_for_atmospheric_window(self, window_name):
        if window_name not in ATMOSPHERIC_WINDOWS:
            raise ValueError(f"Unknown atmospheric window: {window_name}")
        
        window = ATMOSPHERIC_WINDOWS[window_name]
        window_molecules = []

        for molecule in self.molecules.values():
            for molecule_point in molecule.high_and_low_frequencies():
                low_frequency, high_frequency, _ = molecule_point
                for window_low, window_high in window:
                    if window_low <= low_frequency <= window_high or window_low <= high_frequency <= window_high:
                        window_molecules.append(molecule.code)
                        break

        return list(set(window_molecules))

    def filter_for_region_and_intensity(self, region, intensity):
        # Existing implementation
        pass

    def filter_for_region(self, region):
        # Existing implementation
        pass

    def filter_for_intensity(self, intensity, molecule_list=None):
        if molecule_list is None:
            molecule_list = self.molecules.keys()
        
        strong_molecules = []
        for molecule_code in molecule_list:
            molecule = self.molecules[molecule_code]
            for _, _, mol_intensity in molecule.high_and_low_frequencies():
                if mol_intensity >= intensity:
                    strong_molecules.append(molecule_code)
                    break
        return strong_molecules

    # New methods
    def filter_for_wavelength_range(self, start_um, end_um):
        start_cm1 = 10000 / end_um
        end_cm1 = 10000 / start_um
        
        range_molecules = []

        for molecule in self.molecules.values():
            for molecule_point in molecule.high_and_low_frequencies():
                low_frequency, high_frequency, _ = molecule_point
                if start_cm1 <= low_frequency <= end_cm1 or start_cm1 <= high_frequency <= end_cm1:
                    range_molecules.append(molecule.code)
                    break

        return list(set(range_molecules))

    def filter_for_multiple_strong_features(self, threshold, min_features, molecule_list):
        multiple_strong_molecules = []
        for molecule_code in molecule_list:
            molecule = self.molecules[molecule_code]
            strong_features = sum(1 for _, _, intensity in molecule.high_and_low_frequencies() if intensity >= threshold)
            if strong_features >= min_features:
                multiple_strong_molecules.append(molecule_code)
        return multiple_strong_molecules

    def get_strong_features(self, molecule_code, threshold):
        molecule = self.molecules[molecule_code]
        return [(freq, intensity) for low_freq, high_freq, intensity in molecule.high_and_low_frequencies()
                if intensity >= threshold]

