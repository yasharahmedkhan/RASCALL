
import matplotlib
matplotlib.use('TkAgg')

import matplotlib.pyplot as plt

#import seaborn as sns
from . import NIST_spectra
from . import crosssections
import numpy as np

class Plotter:
    colours = ['xkcd:turquoise', 'xkcd:hunter green', 'xkcd:crimson', \
               'xkcd:ochre', 'xkcd:dusty rose', 'xkcd:medium blue', \
               'xkcd:greyish green', 'xkcd:dark peach', 'xkcd:green brown', \
               'xkcd:dark orange', 'xkcd:scarlet', 'xkcd:emerald green', \
               'xkcd:cobalt blue', 'xkcd:neon blue', 'xkcd:evergreen']
    colourIndex = 0

    def average_points(self):
        points = []

        for functional_tuple in self.functionals:
            functional = functional_tuple[0]
            for symmetry in functional.averageSymmetries():
                #                print ('before plotting', symmetry.properties[0].frequency_average())
                for property in symmetry.properties:
                    if property.low != 'UNK':
                        points.append((property.frequency_average(), property.intensity))
                    elif property.low == 'UNK':
                        print ('found unk')

        return points

    def plot_molecule_band_centers(self, molecule):
        # print (len(example.functionals))

        # for functional in molecule
        points = []

        if molecule.isHydrogenMuting:
            functionals =  list(map(lambda functionalTuple: functionalTuple[0], molecule.functionals))
            functionalCodes = list(map(lambda functional: functional.code, functionals))
            print ("found hydrogen muting for molecule ", molecule.code, "with functionals ", functionalCodes)
            self.handleHydrogenMutingMolecule(molecule)
            return

        for functional_tuple in molecule.functionals:
            functional = functional_tuple[0]
            for symmetry in functional.averageSymmetries():
                #                print ('before plotting', symmetry.properties[0].frequency_average())
                for property in symmetry.properties:
                    if property.low != 'UNK':
                        points.append((property.frequency_average(), property.intensity))
                        x = [property.frequency_average()]
                        y = [property.intensity]
                        print (functional.code, symmetry.type, int(property.frequency_average()),\
                            "{0:.2g}".format(property.intensity), self.colours[self.colourIndex])

                        self.setupAppearance(functional, x, y)

                    elif property.low == 'UNK':
                        print ('found unk')

            self.nextColor()
        # plot

        # Plot points
        xs, ys = zip(*points)

        # print ('points', example.points())
        # plt.plot(xs, ys, linestyle='None', marker='o', color='black', linewidth=2)

        # window = range(3250, 3450)
        filtered_list = list(filter(lambda x: 1600 < x[0] < 1900, points))
        # print ('filter:',(filtered_list))

        #print (xs, ys)
       # markerline, stemlines, baseline = plt.stem(xs, ys, '-')
       # plt.setp(baseline, 'color', 'r', 'linewidth', 1)

        # Plot branches
        # for line in example.branches():
        #     x, y = line
        #     markerline, stemlines, baseline = plt.stem(x, y, '-')
        #     plt.setp(baseline, color='r', linewidth=1, marker='None')

        # plt.xlim(1459,1459.5)


        # Plot lines
        # for line in example.lines():
        #    x, y = line

        #    plt.plot(x, y)

    def plot_NIST_spectrum(self, molecule_smile):
        Absorption_Boost = 3
        nu, coef = NIST_spectra.nist_spectrum(molecule_smile)
        plt.plot(nu, coef* Absorption_Boost, label="NIST experimental data")

    def plot_ExoMol_spectrum(self, molecule_smile):
        Absorption_Boost = 3
        hcn_nu = []
        hcn_coef = []
        with open('hcn_exomol.dat', 'rU') as f:
            hcn_exomol = f.readlines()
            for line in hcn_exomol:
                columns = line.strip().split()
                hcn_nu.append(float(columns[0]))
                hcn_coef.append(float(columns[1])*Absorption_Boost)
#        print (hcn_nu, hcn_coef)
        plt.plot(hcn_nu, hcn_coef)

    def plot_ATMOS_crosssections(self, molecule_smile):

        nu, xsec = ATMOS_crosssections.ATMOS_crosssection(molecule_smile)
        plt.plot(nu, xsec, label="ATMOS")



    def show(self, molecule_smile):
        plt.xlabel('Wavenumbers (cm$^{-1}$)', fontsize=16)
        plt.ylabel('Intensity', fontsize=16)
        plt.xlim((0, 4500))
        plt.title(molecule_smile, fontsize=16)

        plt.tick_params(axis='both', which='major', labelsize=12)
        plt.tick_params(axis='both', which='minor', labelsize=12)

        # labels 
        handles, labels = plt.gca().get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        plt.legend(by_label.values(), by_label.keys(), bbox_to_anchor=(1.04, 1), loc="upper left")
                  
        textstr = 'Solid lines are rascall \npredictions and dotted   \nlines are still being\ntested.'

        plt.text(0.78, 0.55, textstr, fontsize=12, transform=plt.gcf().transFigure, verticalalignment="center", horizontalalignment="left", wrap = True)

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
        if reducedCHBondIncidence == 0 :
            shouldMuteCHBondFunctional = True

        for functional_tuple in molecule.functionals:
            functional = functional_tuple[0]

            for symmetry in functional.averageSymmetries():
                #                print ('before plotting', symmetry.properties[0].frequency_average())
                for property in symmetry.properties:
                    if property.low != 'UNK':
                        points.append((property.frequency_average(), property.intensity))
                        x = [property.frequency_average()]
                        if functional.code == chBond and shouldMuteCHBondFunctional:
                            y = [0.001]
                        else:
                            y = [property.intensity]
                        print (functional.code, symmetry.type, int(property.frequency_average()), \
                            "{0:.2g}".format(property.intensity), self.colours[self.colourIndex])

                        self.setupAppearance(functional, x, y)

                    elif property.low == 'UNK':
                        print ('found unk')

            self.nextColor()

    def setupAppearance(self, functional, x, y):
        if functional.source == 'ATMOS':
            markerline, stemlines, baseline = plt.stem(x, y, '--', label= functional.code)
        else:
            markerline, stemlines, baseline = plt.stem(x, y, '-',label= functional.code)

        plt.setp(baseline, visible=False)
        plt.setp(stemlines, 'color', self.colours[self.colourIndex], 'linewidth', 1.5)
        plt.setp(markerline, 'color', self.colours[self.colourIndex], 'linewidth', 1.5)

    def nextColor(self):
        if self.colourIndex == len(self.colours) - 1:
            self.colourIndex = 0
        else:
            self.colourIndex = self.colourIndex + 1

