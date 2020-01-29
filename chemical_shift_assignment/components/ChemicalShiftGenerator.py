import pickle
from timeit import default_timer as timer
from nmr_commons.assignment.ChemicalShiftList import ChemicalShiftList
from nmr_commons.sequence.AminoAcid import AminoAcid
from nmr_commons.sequence.Atom import Atom
from numpy import inf
from nmr_commons.assignment.ChemicalShift import ChemicalShift
from nmr_commons.database.EmpiricalDistribution import EmpiricalDistribution


class ChemicalShiftGenerator:

    def __init__(self, database):
        self.aminoAcidDistributions = []

        if isinstance(database, str):
            self.aminoAcidDistributions = pickle.load(open(database, "rb"))
        else:
            self.generateStatData(database)

    # Save object to file
    def saveAs(self, outputFilePath):
        pickle.dump(self.aminoAcidDistributions, open(outputFilePath, "wb"))

    # Returns empirical distribution
    def getStatisticsByAminoAcid(self, aminoAcidThreeLetterCode):
        return self.aminoAcidDistributions[AminoAcid.find_amino_acid_id(aminoAcidThreeLetterCode)]

    # Visualize empirical distribution of given amino acid
    def visualizeAminoAcid(self, aminoAcidThreeLetterCode):
        id = AminoAcid.find_amino_acid_id(aminoAcidThreeLetterCode)
        for label in AminoAcid.findAminoAcidDefinition(aminoAcidThreeLetterCode):
            print(label)
            self.aminoAcidDistributions[id].plotDistribution(label)

    def visualizeChemicalShift(self, shiftName):

        for code in AminoAcid.THREE_LETTERS_CODES:

            id = AminoAcid.find_amino_acid_id(code)

            if shiftName in AminoAcid.findAminoAcidDefinition(code):
                self.aminoAcidDistributions[id].plotDistribution(shiftName)

    # Function calculates empirical distributions for each amino acid, given BMRB database object
    def generateStatData(self, database):

        # Generate statistics
        for aminoAcid in AminoAcid.THREE_LETTERS_CODES:

            print("Extracting statistics for " + aminoAcid)
            t = timer()

            # Filter amino acids
            result = [elem
                      for record in database.parsed_db
                      for elem in record.sequence
                      if elem.getThreeLettersCode() == aminoAcid]
            print(str(len(result)) + " elements found.")

            # Calculate distribution of each amino acid
            listOfShifts = AminoAcid.findAminoAcidDefinition(aminoAcid)
            listOfShifts = [value for values in listOfShifts.values() for value in values] + list(listOfShifts.keys())
            empiricalDistribution = EmpiricalDistribution(aminoAcid)
            empiricalDistribution.calculateDistribution(result, listOfShifts)

            self.aminoAcidDistributions.append(empiricalDistribution)
            print("Elapsed " + str(timer() - t))

        print("Spectrum generator has been created.")

    # Function generates shift for single amino acid
    def generateShiftsForAminoAcid(self, threeLetterCode, distribution):

        aminoAcidWithShifts = AminoAcid(threeLetterCode)
        listOfShifts = aminoAcidWithShifts.getListOfShifts()
        id = aminoAcidWithShifts.getAminoAcidId()

        listOfShifts = filter(lambda x: Atom.getAtomTypeFromLabel(x) != Atom.ATOM_UNKNOWN, listOfShifts)
        for shift in listOfShifts:
            aminoAcidWithShifts.addChemicalShift(
                ChemicalShift("-", shift, self.aminoAcidDistributions[id].sample(shift, distribution), -1))

        return aminoAcidWithShifts

    # Function generates shifts for amino acid sequence
    def generateShiftsForSequence(self, proteinSequence, distribution):

        outputChemicalShiftlist = ChemicalShiftList()

        # generate chemical shifts for each amino acid
        for aminoAcid in proteinSequence.getSequence():
            generatedShifts = self.generateShiftsForAminoAcid(aminoAcid.getThreeLettersCode(), distribution)
            outputChemicalShiftlist.addAminoAcid(generatedShifts)

        return outputChemicalShiftlist

    def getRangesForAminoAcid(self, aminoAcid):
        ranges = {}
        listOfShifts = filter(lambda x: Atom.getAtomTypeFromLabel(x) != Atom.ATOM_UNKNOWN, aminoAcid.getListOfShifts())

        for shift in listOfShifts:
            ranges[shift] = self.aminoAcidDistributions[aminoAcid.getAminoAcidId()].range(shift)

        return ranges

    def getRangesForSequence(self, proteinSequence):
        ranges = {}
        for aminoAcid in proteinSequence.getSequence():
            partialRanges = self.getRangesForAminoAcid(aminoAcid)
            for key, value in partialRanges.items():
                default = (inf, -inf)
                minRange = min(ranges.get(key, default)[0], value[0])
                maxRange = max(ranges.get(key, default)[1], value[1])
                ranges[key] = (minRange, maxRange)

        return ranges
