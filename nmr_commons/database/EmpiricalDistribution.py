import math
from functools import reduce
import matplotlib.pyplot as plt
import numpy as np
import numpy.random as rnd
from scipy.stats import laplace
from scipy.stats import norm
from nmr_commons.sequence.Atom import Atom


class EmpiricalDistribution:

    LAPLACE = 1
    GAUSS = 2
    UNIFORM = 3
    ALWAYS_MEDIAN = 4
    ALWAYS_MEAN = 5

    # Distribution types
    HISTOGRAM = "hist"
    PDF = "pdf"
    CDF = "cdf"
    BINS = "bins"
    MODEL_PDF = "mpdf"
    MODEL_CDF = "mcdf"

    # Statistics types
    VARIANCE = "var"
    MEAN = "mean"
    MEDIAN = "median"
    BIN_SIZE = "binSize"
    STD = "std"
    BEST = "bEst"

    # Ranges
    N_range = [60, 300]
    C_range = [10, 300]
    CO_range = [10, 300]
    H_range = [4, 12]

    def __init__(self, description=""):
        self.description = description
        self.calculatedDistributions = dict()
        self.numOfBins = 100

    def getLaplaceProbabilityOfChemicalShift(self, shift):
        median = self.getProperty(shift.atomLabel, EmpiricalDistribution.MEDIAN)
        bEst = self.getProperty(shift.atomLabel, EmpiricalDistribution.BEST)
        return laplace.pdf(shift.resonance, median, bEst)


    # Get sample of chemical shift
    def sample(self, shiftName, distribution):
        if distribution == EmpiricalDistribution.UNIFORM:
            ranges = self.getRangeOfChemicalShifts(shiftName)
            return rnd.uniform(*ranges)
        elif distribution == EmpiricalDistribution.LAPLACE:
            median = self.getProperty(shiftName, EmpiricalDistribution.MEDIAN)
            if sum(self.getProperty(shiftName, EmpiricalDistribution.HISTOGRAM)) == 1:
                return median
            else:
                bEst = self.getProperty(shiftName, EmpiricalDistribution.BEST)
                return rnd.laplace(median, bEst)
        elif distribution == EmpiricalDistribution.GAUSS:
            mean = self.getProperty(shiftName, EmpiricalDistribution.MEAN)
            if sum(self.getProperty(shiftName, EmpiricalDistribution.HISTOGRAM)) == 1:
                return mean
            else:
                std = self.getProperty(shiftName, EmpiricalDistribution.STD)
                return rnd.normal(mean, std)
        elif distribution == EmpiricalDistribution.ALWAYS_MEDIAN:
            return self.getProperty(shiftName, EmpiricalDistribution.MEDIAN)
        elif distribution == EmpiricalDistribution.ALWAYS_MEAN:
            return self.getProperty(shiftName, EmpiricalDistribution.MEAN)


    def range(self, shiftName):
        mean = self.getProperty(shiftName, EmpiricalDistribution.MEAN)
        std = self.getProperty(shiftName, EmpiricalDistribution.STD)
        if sum(self.getProperty(shiftName, EmpiricalDistribution.HISTOGRAM)) == 1:
            std = 0
        return mean - 1.5 * std, mean + 1.5 * std

    # Get probability of chemical shift
    def probability(self, shiftName, value):
        median = self.getProperty(shiftName, EmpiricalDistribution.MEDIAN)
        bEst = self.getProperty(shiftName, EmpiricalDistribution.BEST)
        laplace.pdf(value, median, bEst)

    # Caculate empirical distribution
    def calculateDistribution(self, listOfAminoAcids, labels):

        listOfChemicalShifts = [shift
                                for aminoAcid in listOfAminoAcids
                                for shift in aminoAcid.chemicalShifts]

        # For each chemical shift (N, CA, CB, ...)
        for label in filter(lambda x: Atom.getAtomTypeFromLabel(x) != Atom.ATOM_UNKNOWN, labels):

            # Select shifts that match particular label
            listOfSelectedShifts = [float(elem.resonance) for elem in listOfChemicalShifts if elem.atomLabel == label]
            numberOfChemicalShifts = len(listOfSelectedShifts)

            print(str(label) + " " + str(numberOfChemicalShifts))

            # Calculate unnormalized distribution
            histogramRange = self.getRangeOfChemicalShifts(label)
            counts, bins = np.histogram(listOfSelectedShifts, self.numOfBins, histogramRange)

            # Calculate binrange and bin size
            bins = bins[:-1]
            binSize = bins[1] - bins[0]

            # Calculate pdf and cdf
            pdf = list(map(lambda x: x/(numberOfChemicalShifts*binSize), counts))
            cdf = reduce(lambda x, y: x + [x[-1] + y * binSize], pdf, [0])

            # Calculate statistics
            mean = np.mean(listOfSelectedShifts)
            median = np.median(listOfSelectedShifts)
            variance = np.var(listOfSelectedShifts, ddof=1)

            # Laplace distribution
            median = np.median(listOfSelectedShifts)
            bEst = reduce(lambda x, y: x + abs(y - median), listOfSelectedShifts, 0.0)
            bEst = bEst/numberOfChemicalShifts if numberOfChemicalShifts > 0 else None

            # Store all amino acid distributions and statistics in one hashtable
            self.calculatedDistributions[label + "hist"] = counts
            self.calculatedDistributions[label + "bins"] = bins
            self.calculatedDistributions[label + "binSize"] = binSize
            self.calculatedDistributions[label + "pdf"] = pdf
            self.calculatedDistributions[label + "cdf"] = cdf[1:]
            self.calculatedDistributions[label + "mean"] = mean
            self.calculatedDistributions[label + "median"] = median
            self.calculatedDistributions[label + "bEst"] = bEst
            self.calculatedDistributions[label + "var"] = variance
            self.calculatedDistributions[label + "std"] = math.sqrt(variance)
            self.calculatedDistributions[label + "count"] = numberOfChemicalShifts

    # Provide access to properties
    def getProperty(self, shiftLabel, propertyType):
        return self.calculatedDistributions[shiftLabel + propertyType]

    # Define range of empirical distribution
    def getRangeOfChemicalShifts(self, atomType):
        atom = Atom.getAtomTypeFromLabel(atomType)
        if atom == Atom.ATOM_NITROGEN:
            return self.N_range
        elif atom == Atom.ATOM_CO:
            return self.CO_range
        elif atom == Atom.ATOM_CARBON:
            return self.C_range
        elif atom == Atom.ATOM_HYDROGEN:
            return self.H_range

    # Amino acid info plot
    def formatPlot(self, axisObj, title, labelX, labelY):
        titleSize = 11
        labelSize = 10
        ticksSize = 9
        titleSpace = 1.02

        axisObj.set_title(title, fontsize=titleSize, y=titleSpace)
        axisObj.grid()
        axisObj.tick_params(labelsize=ticksSize)
        axisObj.set_xlabel(labelX, fontsize=labelSize)
        axisObj.set_ylabel(labelY, fontsize=labelSize)
        xMin, xMax = axisObj.get_xlim()
        yMin, yMax = axisObj.get_ylim()
        axisObj.set_aspect((xMax - xMin) / (yMax - yMin))

    # Plots empirical distribution
    def plotDistribution(self, label):

        numberOfShifts = self.calculatedDistributions[label + "count"]
        print(self.calculatedDistributions[label + "hist"])
        if numberOfShifts == 0:
            print("Can not visualize " + label + " because number of shifts in database object is zero")
        else:
            print("Visualization based on " + str(numberOfShifts))

            f, axarr = plt.subplots(1, 3)
            f.tight_layout(rect=[0, 0.0, 1, 0.95])
            f.set_facecolor('white')
            plt.suptitle(self.description + " " + label, fontsize=20)

            # xLabels
            xLabels = list(self.calculatedDistributions[label + EmpiricalDistribution.BINS])
            xLabelDelta = xLabels[1] - xLabels[0]

            # Counts
            axarr[0].bar(list(self.calculatedDistributions[label + EmpiricalDistribution.BINS]), list(self.calculatedDistributions[label + EmpiricalDistribution.HISTOGRAM]), xLabelDelta, color="blue")

            # PDF
            axarr[1].bar(list(self.calculatedDistributions[label + EmpiricalDistribution.BINS]),
                            list(self.calculatedDistributions[label + EmpiricalDistribution.PDF]), xLabelDelta,
                            color="blue")
            axarr[1].hold(True)

            mean = self.calculatedDistributions[label + "mean"]
            variance = self.calculatedDistributions[label + "var"]
            median = self.calculatedDistributions[label + "median"]
            bEst = self.calculatedDistributions[label + "bEst"]
            print(str(mean) + " " + str(variance) + " " + str(median) + " " + str(bEst))

            arguments = np.arange(xLabels[0], xLabels[-1], 0.1)

            values = list(map(lambda x: norm.pdf(x, mean, math.sqrt(variance)), arguments))
            axarr[1].plot(arguments, values, linewidth=2, color='r')

            values2 = list(map(lambda x: laplace.pdf(x, median, bEst), arguments))
            axarr[1].plot(arguments, values2, linewidth=2, color='g')

            # CDF
            axarr[2].bar(list(self.calculatedDistributions[label + EmpiricalDistribution.BINS]),
                            list(self.calculatedDistributions[label + EmpiricalDistribution.CDF]), xLabelDelta,
                            color="blue")
            axarr[2].hold(True)
            arguments = np.arange(xLabels[0], xLabels[-1], 0.1)

            values = list(map(lambda x: norm.cdf(x, mean, math.sqrt(variance)), arguments))
            axarr[2].plot(arguments, values, linewidth=2, color='r')

            values2 = list(map(lambda x: laplace.cdf(x, median, bEst), arguments))
            axarr[2].plot(arguments, values2, linewidth=2, color='g')

            # Format plots
            self.formatPlot(axarr[0], "Histogram", "ppm", "counts")
            self.formatPlot(axarr[1], "PDF", "ppm", "prob")
            self.formatPlot(axarr[2], "CDF", "ppm", "prob")
            plt.draw()
