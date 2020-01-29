from math import pi
import numpy as np
from nmr_commons.spectrum.Spectrum import Spectrum
from nmr_commons.peak_picking.PeakList import PeakList


class PeakListGenParams:
    def __init__(self):
        self.extraTrueProb, self.extraTrueMax = 0, 0
        self.noisyProb, self.noisyMax = 0, 0
        self.incompleteProb = 0
        self.noisyRndGenParam = 1 / 65.
        self.noisyRndGen = lambda: 2 / pi * np.arctan(np.random.exponential(scale=self.noisyRndGenParam))
        self.trueRndGen = lambda: 1 - self.noisyRndGen()

    def copy(self):
        cp = PeakListGenParams().extra_true(self.extraTrueMax, self.extraTrueProb)\
            .incomplete(self.incompleteProb).noisy(self.noisyMax, self.noisyProb)
        cp.noisyRndGenParam = self.noisyRndGenParam
        cp.noisyRndGen = self.noisyRndGen
        cp.trueRndGen = self.trueRndGen
        return cp

    def extra_true(self, max_extra_per_one_real, prob=1):
        self.extraTrueProb, self.extraTrueMax = prob, max_extra_per_one_real
        return self

    def noisy(self, max_noisy_per_one_real, prob=1):
        self.noisyProb, self.noisyMax = prob, max_noisy_per_one_real
        return self

    def incomplete(self, prob):
        self.incompleteProb = prob
        return self


# CLASS FUNCTIONS AND PROPERTIES:
#   - Generates peak positions, given chemical shift assignment
class PeakListGenerator:

    @staticmethod
    def generatePeak(chemShiftList, peakDefinition, aminoAcidId):
        peakCoordinates = []
        containsAllShifts = True
        missingShifts = ""
        for shiftModifier, shiftName in peakDefinition:

            if aminoAcidId + shiftModifier < 0:
                missingShifts += shiftName + " "
                continue

            aminoAcidObj = chemShiftList.getAminoAcidById(aminoAcidId + shiftModifier)

            containsAllShifts = containsAllShifts and aminoAcidObj.containsShift(shiftName)
            shift = aminoAcidObj.getChemicalShiftByLabel(shiftName)

            if shift is not None:
                peakCoordinates.append(shift.resonance)
            else:
                missingShifts += shiftName + " "

        return peakCoordinates, containsAllShifts, missingShifts

    #
    # Generates list of peaks given protein sequence and spectrum definition
    @staticmethod
    def generate_peak_list(mainChemShiftList, spectrumDefinition, params=None, extraChemShiftLists=[], verbose=True):
        if params is None:
            params = PeakListGenParams()
        if params.extraTrueMax > len(extraChemShiftLists):
            raise AttributeError("ExtraChemShiftLists size is smaller than params.extraTrueMax.")

        outputPeakList = PeakList([key['label'] for key in spectrumDefinition['axes']], items=[])
        numberOfIncompletePeaks = 0
        numberOfCompletePeaks = 0

        # iterate over all amino acids
        for aminoAcidId in range(mainChemShiftList.get_number_of_amino_acids()):
            peakDefs = spectrumDefinition['peakDefinition'] if 'peakDefinition' in spectrumDefinition \
                else Spectrum.get_peak_definition(mainChemShiftList.getAminoAcidById(aminoAcidId - 1),
                                                  mainChemShiftList.getAminoAcidById(aminoAcidId), spectrumDefinition)

            for peakDefinition in peakDefs:
                peak, defContainsShifts, missingShifts = PeakListGenerator.generatePeak(mainChemShiftList,
                                                                                        peakDefinition, aminoAcidId)

                # if peak generation was successful
                if len(peak) == len(peakDefinition):
                    numberOfCompletePeaks += 1
                    if verbose: print("--- New peak addition phase...")
                    if np.random.uniform() > params.incompleteProb:
                        outputPeakList.add_peak(peak, params.trueRndGen())
                        if verbose: print("Classic peak added. " + str(outputPeakList.get_list()[-1]))
                    else:
                        if verbose: print("Classic peak rejected. " + str(peak))

                    for i in range(params.extraTrueMax):
                        if np.random.uniform() < params.extraTrueProb:
                            extraPeak, _, _ = PeakListGenerator.generatePeak(extraChemShiftLists[i],
                                                                             peakDefinition, aminoAcidId)
                            outputPeakList.add_peak(extraPeak, params.trueRndGen())
                            if verbose: print("Extra peak added. " + str(outputPeakList.get_list()[-1]))

                    for i in range(params.noisyMax):
                        if np.random.uniform() < params.noisyProb:
                            noisyPeak = [np.random.uniform(axis['PPMRange'][0], axis['PPMRange'][1])
                                         for axis in spectrumDefinition['axes']]
                            outputPeakList.add_peak(noisyPeak, params.noisyRndGen())
                            if verbose: print("Noisy peak added. " + str(outputPeakList.get_list()[-1]))

                elif defContainsShifts:
                    numberOfIncompletePeaks += 1
                    if verbose:
                        print("! Incomplete assignment for: ")
                        print("   " + mainChemShiftList.getAminoAcidById(aminoAcidId).toString(full=True))
                        print("   Missing shifts: " + missingShifts)

        return outputPeakList, numberOfIncompletePeaks, numberOfCompletePeaks
