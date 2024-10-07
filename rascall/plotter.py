import matplotlib
matplotlib.use('TkAgg')

import matplotlib.pyplot as plt
import numpy as np
from . import NIST_spectra
from . import crosssections

class Plotter:
    colours = ['xkcd:turquoise', 'xkcd:hunter green', 'xkcd:crimson',
               'xkcd:ochre', 'xkcd:dusty rose', 'xkcd:medium blue',
               'xkcd:greyish green', 'xkcd:dark peach', 'xkcd:green brown',
               'xkcd:dark orange', 'xkcd:scarlet', 'xkcd:emerald green',
               'xkcd:cobalt blue', 'xkcd:neon blue', 'xkcd:evergreen']
    colourIndex = 0

    def get_molecule_band_centers(self, molecule):
        points = {}
        
        if molecule.isHydrogenMuting:
            functionals = list(map(lambda functionalTuple: functionalTuple[0], molecule.functionals))
            functionalCodes = list(map(lambda functional: functional.code, functionals))
            print("found hydrogen muting for molecule ", molecule.code, "with functionals ", functionalCodes)
            return self.get_hydrogen_muting_molecule_data(molecule)

        for functional_tuple in molecule.functionals:
            functional = functional_tuple[0]
            points[functional.code] = []
            for symmetry in functional.averageSymmetries():
                for property in symmetry.properties:
                    if property.low != 'UNK':
                        points[functional.code].append((property.frequency_average(), property.intensity))
                        print(functional.code, symmetry.type, int(property.frequency_average()),
                              "{0:.2g}".format(property.intensity), self.colours[self.colourIndex])
                    elif property.low == 'UNK':
                        print('found unk')

            self.nextColor()

        return points

    def get_hydrogen_muting_molecule_data(self, molecule):
        points = {}
        tripleBond = '[H]C#[!#1]'
        chBond = 'C[H]'
        tripleBondIncidence = 0
        chBondIncidence = 0
        shouldMuteCHBondFunctional = False

        for functionalTuple in molecule.functionals:
            functional = functionalTuple[0]
            functionalIncidence = functionalTuple[1]

            if functional.code == tripleBond:
                tripleBondIncidence = functionalIncidence
            elif functional.code == chBond:
                chBondIncidence = functionalIncidence

        reducedCHBondIncidence = max(0, tripleBondIncidence - chBondIncidence)
        if reducedCHBondIncidence == 0:
            shouldMuteCHBondFunctional = True

        for functional_tuple in molecule.functionals:
            functional = functional_tuple[0]
            points[functional.code] = []

            for symmetry in functional.averageSymmetries():
                for property in symmetry.properties:
                    if property.low != 'UNK':
                        intensity = property.intensity if not (functional.code == chBond and shouldMuteCHBondFunctional) else 0.001
                        points[functional.code].append((property.frequency_average(), intensity))
                        print(functional.code, symmetry.type, int(property.frequency_average()),
                              "{0:.2g}".format(intensity), self.colours[self.colourIndex])
                    elif property.low == 'UNK':
                        print('found unk')

            self.nextColor()

        return points

    def get_NIST_spectrum(self, molecule_smile):
        Absorption_Boost = 3
        nu, coef = NIST_spectra.nist_spectrum(molecule_smile)
        return list(zip(nu, coef * Absorption_Boost))

    def nextColor(self):
        if self.colourIndex == len(self.colours) - 1:
            self.colourIndex = 0
        else:
            self.colourIndex = self.colourIndex + 1

    # Existing methods for terminal-based plotting
    def plot_molecule_band_centers(self, molecule):
        points = []

        if molecule.isHydrogenMuting:
            functionals = list(map(lambda functionalTuple: functionalTuple[0], molecule.functionals))
            functionalCodes = list(map(lambda functional: functional.code, functionals))
            print("found hydrogen muting for molecule ", molecule.code, "with functionals ", functionalCodes)
            self.handleHydrogenMutingMolecule(molecule)
            return

        for functional_tuple in molecule.functionals:
            functional = functional_tuple[0]
            for symmetry in functional.averageSymmetries():
                for property in symmetry.properties:
                    if property.low != 'UNK':
                        points.append((property.frequency_average(), property.intensity))
                        x = [property.frequency_average()]
                        y = [property.intensity]
                        print(functional.code, symmetry.type, int(property.frequency_average()),
                            "{0:.2g}".format(property.intensity), self.colours[self.colourIndex])

                        self.setupAppearance(functional, x, y)

                    elif property.low == 'UNK':
                        print('found unk')

            self.nextColor()

        xs, ys = zip(*points)

    def plot_NIST_spectrum(self, molecule_smile):
        Absorption_Boost = 3
        nu, coef = NIST_spectra.nist_spectrum(molecule_smile)
        plt.plot(nu, coef * Absorption_Boost, label="NIST experimental data")

    def show(self, molecule_smile):
        plt.xlabel('Wavenumbers (cm$^{-1}$)', fontsize=16)
        plt.ylabel('Intensity', fontsize=16)
        plt.xlim((0, 4500))
        plt.title(molecule_smile, fontsize=16)

        plt.tick_params(axis='both', which='major', labelsize=12)
        plt.tick_params(axis='both', which='minor', labelsize=12)

        handles, labels = plt.gca().get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        plt.legend(by_label.values(), by_label.keys(), bbox_to_anchor=(1.04, 1), loc="upper left")
                  
        textstr = 'Solid lines are rascall \npredictions and dotted   \nlines are still being\ntested.'

        plt.text(0.78, 0.55, textstr, fontsize=12, transform=plt.gcf().transFigure, verticalalignment="center", horizontalalignment="left", wrap=True)

        fig = matplotlib.pyplot.gcf()
        fig.set_size_inches(10, 6)
        plt.tight_layout()       
        plt.show()

    def handleHydrogenMutingMolecule(self, molecule):
        points = []
        tripleBond = '[H]C#[!#1]'
        chBond = 'C[H]'
        tripleBondIncidence = 0
        chBondIncidence = 0
        shouldMuteCHBondFunctional = False

        for functionalTuple in molecule.functionals:
            functional = functionalTuple[0]
            functionalIncidence = functionalTuple[1]

            if functional.code == tripleBond:
                tripleBondIncidence = functionalIncidence
            elif functional.code == chBond:
                chBondIncidence = functionalIncidence

        reducedCHBondIncidence = max(0, tripleBondIncidence - chBondIncidence)
        if reducedCHBondIncidence == 0:
            shouldMuteCHBondFunctional = True

        for functional_tuple in molecule.functionals:
            functional = functional_tuple[0]

            for symmetry in functional.averageSymmetries():
                for property in symmetry.properties:
                    if property.low != 'UNK':
                        points.append((property.frequency_average(), property.intensity))
                        x = [property.frequency_average()]
                        if functional.code == chBond and shouldMuteCHBondFunctional:
                            y = [0.001]
                        else:
                            y = [property.intensity]
                        print(functional.code, symmetry.type, int(property.frequency_average()),
                            "{0:.2g}".format(property.intensity), self.colours[self.colourIndex])

                        self.setupAppearance(functional, x, y)

                    elif property.low == 'UNK':
                        print('found unk')

            self.nextColor()

    def setupAppearance(self, functional, x, y):
        if functional.source == 'ATMOS':
            markerline, stemlines, baseline = plt.stem(x, y, '--', label=functional.code)
        else:
            markerline, stemlines, baseline = plt.stem(x, y, '-', label=functional.code)

        plt.setp(baseline, visible=False)
        plt.setp(stemlines, 'color', self.colours[self.colourIndex], 'linewidth', 1.5)
        plt.setp(markerline, 'color', self.colours[self.colourIndex], 'linewidth', 1.5)




