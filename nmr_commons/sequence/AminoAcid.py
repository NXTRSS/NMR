import warnings
from collections import OrderedDict
from nmr_commons.sequence.Atom import Atom
from nmr_commons.assignment.ChemicalShift import ChemicalShift


class AminoAcid:

    # AMINO ACID DEFINITIONS
    # ARTIFICIAL_AMINOACID = OrderedDict([
    #     (Atom.CARBON, []),
    #     (Atom.NITROGEN, [Atom.HYDROGEN]),
    #     (Atom.CARBON_ALPHA, [Atom.HYDROGEN_ALPHA, Atom.HYDROGEN_ALPHA2, Atom.HYDROGEN_ALPHA3]),
    #     (Atom.CARBON_BETA, [Atom.HYDROGEN_BETA, Atom.HYDROGEN_BETA2, Atom.HYDROGEN_BETA3]),
    #     (Atom.CARBON_GAMMA, [Atom.HYDROGEN_GAMMA, Atom.HYDROGEN_GAMMA2, Atom.HYDROGEN_GAMMA3]),
    #     (Atom.CARBON_GAMMA1, [Atom.HYDROGEN_GAMMA1, Atom.HYDROGEN_GAMMA12, Atom.HYDROGEN_GAMMA13]),
    #     (Atom.CARBON_GAMMA2, [Atom.HYDROGEN_GAMMA2]),
    #     (Atom.OXYGEN_GAMMA, [Atom.HYDROGEN_GAMMA]),
    #     (Atom.OXYGEN_GAMMA1, [Atom.HYDROGEN_GAMMA1]),
    #     (Atom.CARBON_DELTA, [Atom.HYDROGEN_DELTA2, Atom.HYDROGEN_DELTA3]),
    #     (Atom.CARBON_DELTA1, [Atom.HYDROGEN_DELTA1]),
    #     (Atom.CARBON_DELTA2, [Atom.HYDROGEN_DELTA2]),
    #     (Atom.NITROGEN_DELTA1, [Atom.HYDROGEN_DELTA1]),
    #     (Atom.NITROGEN_DELTA2, [Atom.HYDROGEN_DELTA21, Atom.HYDROGEN_DELTA22]),
    #     (Atom.OXYGEN_DELTA, [Atom.HYDROGEN_DELTA2]),
    #     (Atom.CARBON_EPSILON, [Atom.HYDROGEN_EPSILON, Atom.HYDROGEN_EPSILON2, Atom.HYDROGEN_EPSILON3]),
    #     (Atom.CARBON_EPSILON1, [Atom.HYDROGEN_EPSILON1]),
    #     (Atom.CARBON_EPSILON2, [Atom.HYDROGEN_EPSILON2]),
    #     (Atom.CARBON_EPSILON3, [Atom.HYDROGEN_EPSILON3]),
    #     (Atom.NITROGEN_EPSILON, [Atom.HYDROGEN_EPSILON]),
    #     (Atom.NITROGEN_EPSILON1, [Atom.HYDROGEN_EPSILON1]),
    #     (Atom.NITROGEN_EPSILON2, [Atom.HYDROGEN_EPSILON2, Atom.HYDROGEN_EPSILON21, Atom.HYDROGEN_EPSILON22]),
    #     (Atom.OXYGEN_EPSILON, [Atom.HYDROGEN_EPSILON2]),
    #     (Atom.NITROGEN_ZETA, [Atom.HYDROGEN_ZETA]),
    #     (Atom.CARBON_ZETA, [Atom.HYDROGEN_ZETA]),
    #     (Atom.CARBON_ZETA2, [Atom.HYDROGEN_ZETA2]),
    #     (Atom.CARBON_ZETA3, [Atom.HYDROGEN_ZETA3]),
    #     (Atom.NITROGEN_NH1, [Atom.HYDROGEN_H11, Atom.HYDROGEN_H12]),
    #     (Atom.NITROGEN_NH2, [Atom.HYDROGEN_H21, Atom.HYDROGEN_H22]),
    #     (Atom.CARBON_H2, [Atom.HYDROGEN_HH2]),
    #     (Atom.OXYGEN_ETA, [Atom.HYDROGEN_HH]),
    # ])

    ALANINE = OrderedDict([
        (Atom.CARBON, []),
        (Atom.NITROGEN, [Atom.HYDROGEN]),
        (Atom.CARBON_ALPHA, [Atom.HYDROGEN_ALPHA]),
        (Atom.CARBON_BETA, [Atom.HYDROGEN_BETA])
    ])

    ARGININE = OrderedDict([
        (Atom.CARBON, []),
        (Atom.NITROGEN, [Atom.HYDROGEN]),
        (Atom.CARBON_ALPHA, [Atom.HYDROGEN_ALPHA]),
        (Atom.CARBON_BETA, [Atom.HYDROGEN_BETA2, Atom.HYDROGEN_BETA3]),
        (Atom.CARBON_GAMMA, [Atom.HYDROGEN_GAMMA2, Atom.HYDROGEN_GAMMA3]),
        (Atom.CARBON_DELTA, [Atom.HYDROGEN_DELTA2, Atom.HYDROGEN_DELTA3]),
        (Atom.NITROGEN_EPSILON, [Atom.HYDROGEN_EPSILON]),
        (Atom.CARBON_ZETA, []),
        (Atom.NITROGEN_NH1, [Atom.HYDROGEN_H11, Atom.HYDROGEN_H12]),
        (Atom.NITROGEN_NH2, [Atom.HYDROGEN_H21, Atom.HYDROGEN_H22])
    ])

    ASPARTIC_ACID = OrderedDict([
        (Atom.CARBON, []),
        (Atom.NITROGEN, [Atom.HYDROGEN]),
        (Atom.CARBON_ALPHA, [Atom.HYDROGEN_ALPHA]),
        (Atom.CARBON_BETA, [Atom.HYDROGEN_BETA2, Atom.HYDROGEN_BETA3]),
        (Atom.CARBON_GAMMA, []),
        (Atom.OXYGEN_DELTA, [Atom.HYDROGEN_DELTA2])
    ])

    ASPARAGINE = OrderedDict([
        (Atom.CARBON, []),
        (Atom.NITROGEN, [Atom.HYDROGEN]),
        (Atom.CARBON_ALPHA, [Atom.HYDROGEN_ALPHA]),
        (Atom.CARBON_BETA, [Atom.HYDROGEN_BETA2, Atom.HYDROGEN_BETA3]),
        (Atom.CARBON_GAMMA, []),
        (Atom.NITROGEN_DELTA2, [Atom.HYDROGEN_DELTA21, Atom.HYDROGEN_DELTA22])
    ])

    CYSTEINE = OrderedDict([
        (Atom.CARBON, []),
        (Atom.NITROGEN, [Atom.HYDROGEN]),
        (Atom.CARBON_ALPHA, [Atom.HYDROGEN_ALPHA]),
        (Atom.CARBON_BETA, [Atom.HYDROGEN_BETA2, Atom.HYDROGEN_BETA3]),
        (Atom.SULFUR_GAMMA, [Atom.HYDROGEN_GAMMA])
    ])

    GLUTAMIC_ACID = OrderedDict([
        (Atom.CARBON, []),
        (Atom.NITROGEN, [Atom.HYDROGEN]),
        (Atom.CARBON_ALPHA, [Atom.HYDROGEN_ALPHA]),
        (Atom.CARBON_BETA, [Atom.HYDROGEN_BETA2, Atom.HYDROGEN_BETA3]),
        (Atom.CARBON_GAMMA, [Atom.HYDROGEN_GAMMA2, Atom.HYDROGEN_GAMMA3]),
        (Atom.CARBON_DELTA, []),
        (Atom.OXYGEN_EPSILON, [Atom.HYDROGEN_EPSILON2])
    ])

    GLUTAMINE = OrderedDict([
        (Atom.CARBON, []),
        (Atom.NITROGEN, [Atom.HYDROGEN]),
        (Atom.CARBON_ALPHA, [Atom.HYDROGEN_ALPHA]),
        (Atom.CARBON_BETA, [Atom.HYDROGEN_BETA2, Atom.HYDROGEN_BETA3]),
        (Atom.CARBON_GAMMA, [Atom.HYDROGEN_GAMMA2, Atom.HYDROGEN_GAMMA3]),
        (Atom.CARBON_DELTA, []),
        (Atom.NITROGEN_EPSILON2, [Atom.HYDROGEN_EPSILON21, Atom.HYDROGEN_EPSILON22])
    ])

    GLYCINE = OrderedDict([
        (Atom.CARBON, []),
        (Atom.NITROGEN, [Atom.HYDROGEN]),
        (Atom.CARBON_ALPHA, [Atom.HYDROGEN_ALPHA2, Atom.HYDROGEN_ALPHA3])
    ])

    HISTIDINE = OrderedDict([
        (Atom.CARBON, []),
        (Atom.NITROGEN, [Atom.HYDROGEN]),
        (Atom.CARBON_ALPHA, [Atom.HYDROGEN_ALPHA]),
        (Atom.CARBON_BETA, [Atom.HYDROGEN_BETA2, Atom.HYDROGEN_BETA3]),
        (Atom.CARBON_GAMMA, []),
        (Atom.NITROGEN_DELTA1, [Atom.HYDROGEN_DELTA1]),
        (Atom.CARBON_DELTA2, [Atom.HYDROGEN_DELTA2]),
        (Atom.CARBON_EPSILON1, [Atom.HYDROGEN_EPSILON1]),
        (Atom.NITROGEN_EPSILON2, [Atom.HYDROGEN_EPSILON2])
    ])

    ISOLEUCINE = OrderedDict([
        (Atom.CARBON, []),
        (Atom.NITROGEN, [Atom.HYDROGEN]),
        (Atom.CARBON_ALPHA, [Atom.HYDROGEN_ALPHA]),
        (Atom.CARBON_BETA, [Atom.HYDROGEN_BETA]),
        (Atom.CARBON_GAMMA1, [Atom.HYDROGEN_GAMMA12, Atom.HYDROGEN_GAMMA13]),
        (Atom.CARBON_GAMMA2, [Atom.HYDROGEN_GAMMA2]),
        (Atom.CARBON_DELTA1, [Atom.HYDROGEN_DELTA1])
    ])

    LEUCINE = OrderedDict([
        (Atom.CARBON, []),
        (Atom.NITROGEN, [Atom.HYDROGEN]),
        (Atom.CARBON_ALPHA, [Atom.HYDROGEN_ALPHA]),
        (Atom.CARBON_BETA, [Atom.HYDROGEN_BETA2, Atom.HYDROGEN_BETA3]),
        (Atom.CARBON_GAMMA, [Atom.HYDROGEN_GAMMA]),
        (Atom.CARBON_DELTA1, [Atom.HYDROGEN_DELTA1]),
        (Atom.CARBON_DELTA2, [Atom.HYDROGEN_DELTA2])
    ])

    LYSINE = OrderedDict([
        (Atom.CARBON, []),
        (Atom.NITROGEN, [Atom.HYDROGEN]),
        (Atom.CARBON_ALPHA, [Atom.HYDROGEN_ALPHA]),
        (Atom.CARBON_BETA, [Atom.HYDROGEN_BETA2, Atom.HYDROGEN_BETA3]),
        (Atom.CARBON_GAMMA, [Atom.HYDROGEN_GAMMA2, Atom.HYDROGEN_GAMMA3]),
        (Atom.CARBON_DELTA, [Atom.HYDROGEN_DELTA2, Atom.HYDROGEN_DELTA3]),
        (Atom.CARBON_EPSILON, [Atom.HYDROGEN_EPSILON2, Atom.HYDROGEN_EPSILON3]),
        (Atom.NITROGEN_ZETA, [Atom.HYDROGEN_ZETA])
    ])

    METHIONINE = OrderedDict([
        (Atom.CARBON, []),
        (Atom.NITROGEN, [Atom.HYDROGEN]),
        (Atom.CARBON_ALPHA, [Atom.HYDROGEN_ALPHA]),
        (Atom.CARBON_BETA, [Atom.HYDROGEN_BETA2, Atom.HYDROGEN_BETA3]),
        (Atom.CARBON_GAMMA, [Atom.HYDROGEN_GAMMA2, Atom.HYDROGEN_GAMMA3]),
        (Atom.SULFUR_DELTA, []),
        (Atom.CARBON_EPSILON, [Atom.HYDROGEN_EPSILON])
    ])

    PHENYLALANINE = OrderedDict([
        (Atom.CARBON, []),
        (Atom.NITROGEN, [Atom.HYDROGEN]),
        (Atom.CARBON_ALPHA, [Atom.HYDROGEN_ALPHA]),
        (Atom.CARBON_BETA, [Atom.HYDROGEN_BETA2, Atom.HYDROGEN_BETA3]),
        (Atom.CARBON_GAMMA, []),
        (Atom.CARBON_DELTA1, [Atom.HYDROGEN_DELTA1]),
        (Atom.CARBON_DELTA2, [Atom.HYDROGEN_DELTA2]),
        (Atom.CARBON_EPSILON1, [Atom.HYDROGEN_EPSILON1]),
        (Atom.CARBON_EPSILON2, [Atom.HYDROGEN_EPSILON2]),
        (Atom.CARBON_ZETA, [Atom.HYDROGEN_ZETA])
    ])

    PROLINE = OrderedDict([
        (Atom.CARBON, []),
        (Atom.CARBON_ALPHA, [Atom.HYDROGEN_ALPHA]),
        (Atom.CARBON_BETA, [Atom.HYDROGEN_BETA2, Atom.HYDROGEN_BETA3]),
        (Atom.CARBON_GAMMA, [Atom.HYDROGEN_GAMMA2, Atom.HYDROGEN_GAMMA3]),
        (Atom.CARBON_DELTA, [Atom.HYDROGEN_DELTA2, Atom.HYDROGEN_DELTA3]),
        (Atom.NITROGEN, [])
    ])

    SERINE = OrderedDict([
        (Atom.CARBON, []),
        (Atom.NITROGEN, [Atom.HYDROGEN]),
        (Atom.CARBON_ALPHA, [Atom.HYDROGEN_ALPHA]),
        (Atom.CARBON_BETA, [Atom.HYDROGEN_BETA2, Atom.HYDROGEN_BETA3]),
        (Atom.OXYGEN_GAMMA, [Atom.HYDROGEN_GAMMA])
    ])

    THREONINE = OrderedDict([
        (Atom.CARBON, []),
        (Atom.NITROGEN, [Atom.HYDROGEN]),
        (Atom.CARBON_ALPHA, [Atom.HYDROGEN_ALPHA]),
        (Atom.CARBON_BETA, [Atom.HYDROGEN_BETA]),
        (Atom.OXYGEN_GAMMA1, [Atom.HYDROGEN_GAMMA1]),
        (Atom.CARBON_GAMMA2, [Atom.HYDROGEN_GAMMA2])
    ])

    TRYPTOPHAN = OrderedDict([
        (Atom.CARBON, []),
        (Atom.NITROGEN, [Atom.HYDROGEN]),
        (Atom.CARBON_ALPHA, [Atom.HYDROGEN_ALPHA]),
        (Atom.CARBON_BETA, [Atom.HYDROGEN_BETA2, Atom.HYDROGEN_BETA3]),
        (Atom.CARBON_GAMMA, []),
        (Atom.CARBON_DELTA1, [Atom.HYDROGEN_DELTA1]),
        (Atom.CARBON_DELTA2, []),
        (Atom.NITROGEN_EPSILON1, [Atom.HYDROGEN_EPSILON1]),
        (Atom.CARBON_EPSILON2, []),
        (Atom.CARBON_EPSILON3, [Atom.HYDROGEN_EPSILON3]),
        (Atom.CARBON_ZETA2, [Atom.HYDROGEN_ZETA2]),
        (Atom.CARBON_ZETA3, [Atom.HYDROGEN_ZETA3]),
        (Atom.CARBON_H2, [Atom.HYDROGEN_HH2])
    ])

    TYROSINE = OrderedDict([
        (Atom.CARBON, []),
        (Atom.NITROGEN, [Atom.HYDROGEN]),
        (Atom.CARBON_ALPHA, [Atom.HYDROGEN_ALPHA]),
        (Atom.CARBON_BETA, [Atom.HYDROGEN_BETA2, Atom.HYDROGEN_BETA3]),
        (Atom.CARBON_GAMMA, []),
        (Atom.CARBON_DELTA1, [Atom.HYDROGEN_DELTA1]),
        (Atom.CARBON_DELTA2, [Atom.HYDROGEN_DELTA2]),
        (Atom.CARBON_EPSILON1, [Atom.HYDROGEN_EPSILON1]),
        (Atom.CARBON_EPSILON2, [Atom.HYDROGEN_EPSILON2]),
        (Atom.CARBON_ZETA, []),
        (Atom.OXYGEN_ETA, [Atom.HYDROGEN_HH])
    ])

    VALINE = OrderedDict([
        (Atom.CARBON, []),
        (Atom.NITROGEN, [Atom.HYDROGEN]),
        (Atom.CARBON_ALPHA, [Atom.HYDROGEN_ALPHA]),
        (Atom.CARBON_BETA, [Atom.HYDROGEN_BETA]),
        (Atom.CARBON_GAMMA1, [Atom.HYDROGEN_GAMMA1]),
        (Atom.CARBON_GAMMA2, [Atom.HYDROGEN_GAMMA2])
    ])

    UNKNOWN = OrderedDict([])

    AMINO_ACID_DEFINITIONS = [
        ALANINE, ARGININE, ASPARAGINE, ASPARTIC_ACID, CYSTEINE,
        GLUTAMIC_ACID, GLUTAMINE, GLYCINE, HISTIDINE, ISOLEUCINE,
        LEUCINE, LYSINE, METHIONINE, PHENYLALANINE, PROLINE, SERINE, THREONINE,
        TRYPTOPHAN, TYROSINE, VALINE, UNKNOWN]  # , ARTIFICIAL_AMINOACID]

    # AMINO ACID NAMES
    THREE_LETTERS_CODES = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE", "LEU", "LYS", "MET",
                           "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "UUU"]
    ONE_LETTER_CODES = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y",
                        "V", "U"]
    NAMES = [
        'alanine', 'arginine', 'asparagine', 'aspartic_acid', 'cysteine',
        'glutamic_acid', 'glutamine', 'glycine', 'histidine', 'isoleucine',
        'leucine', 'lysine', 'methionine', 'phenylalanine', 'proline', 'serine', 'threonine',
        'tryptophan', 'tyrosine', 'valine', 'unknown'
    ]

    # WEIGHTS = [
    #     89.1, 174.2, 132.1, 133.1, 121.2, 147.1, 146.2, 75.1, 155.2, 131.2, 131.2, 146.2, 149.2, 165.2, 115.1, 105.1, 119.1, 204.2, 181.2, 117.1, 0
    # ]

    WEIGHTS = [
        71.037, 156.18, 114.10, 115.08, 103.14, 129.11,128.058, 57.0513, 137.1393, 113.1576,
        113.1576, 128.1723, 131.1961, 147.1739, 97.1152, 87.0773, 101.1039,
        186.2099, 163.1733, 99.1311, 0
    ]

    NUMBER_OF_AMINO_ACIDS = len(THREE_LETTERS_CODES)

    @staticmethod
    def getAllShiftNamesFromAminoAcidDefinition(aminoAcidDefinition):

        list = []

        for key in aminoAcidDefinition.items():
            list.append(key[0])
            list = list + key[1]

        return list

    def getAllShiftNamesFromDefinition(self):
        return self.getAllShiftNamesFromAminoAcidDefinition(self.definition)

    @staticmethod
    def find_amino_acid_id(amino_acid_identifier):

        # For one letters code
        if len(amino_acid_identifier) == 1:
            if amino_acid_identifier not in AminoAcid.ONE_LETTER_CODES:
                amino_acid_identifier = "U"
            return AminoAcid.ONE_LETTER_CODES.index(amino_acid_identifier)

        # For three letters code
        elif len(amino_acid_identifier) == 3:
            if amino_acid_identifier not in AminoAcid.THREE_LETTERS_CODES:
                amino_acid_identifier = "UUU"
            return AminoAcid.THREE_LETTERS_CODES.index(amino_acid_identifier)

        # Unknown
        else:
            raise IndexError("Can not interpret " + str(amino_acid_identifier) + " as amino acid identifier.")

    @staticmethod
    def findAminoAcidDefinition(aminoAcidIdentifier):
        if type(aminoAcidIdentifier) == int:
            return AminoAcid.AMINO_ACID_DEFINITIONS[aminoAcidIdentifier]
        else:
            return AminoAcid.AMINO_ACID_DEFINITIONS[AminoAcid.find_amino_acid_id(aminoAcidIdentifier)]

    def __init__(self, aminoAcidIdentifier):
        self.chemicalShifts = []
        self.aminoAcidId = AminoAcid.find_amino_acid_id(aminoAcidIdentifier)
        self.definition = AminoAcid.findAminoAcidDefinition(self.getAminoAcidId())
        self.defined_shift_list = [value for values in self.definition.values() for value in values] + list(self.definition.keys())

    def get_weight_in_Da(self):
        return self.WEIGHTS[self.aminoAcidId]

    def serialize_z(self, listOfConsideredShifts):
        listOfShifts = AminoAcid.getAllShiftNamesFromAminoAcidDefinition(self.definition)

        list = []
        for shift in listOfShifts:

            if Atom.getAtomTypeFromLabel(shift) != Atom.ATOM_UNKNOWN and (shift in listOfConsideredShifts):
                ch = self.getChemicalShiftByLabel(shift)
                list.append( ch.resonance if not ch == None else None )

        return list

    def getShiftOrder(self):
        listOfShifts = AminoAcid.getAllShiftNamesFromAminoAcidDefinition(self.definition)
        return [item for item in listOfShifts if Atom.getAtomTypeFromLabel(item) != Atom.ATOM_UNKNOWN]

    def getOneLetterCode(self):
        return self.ONE_LETTER_CODES[self.aminoAcidId]

    def getThreeLettersCode(self):
        return self.THREE_LETTERS_CODES[self.aminoAcidId]

    def getAminoAcidId(self):
        return self.aminoAcidId

    def getChemicalShifts(self):
        return self.chemicalShifts

    def getChemicalShiftByLabel(self, label):
        for chemicalShift in self.chemicalShifts:
            if chemicalShift.atomLabel == label:
                return chemicalShift

    def addChemicalShift(self, shift):

        if isinstance(shift, ChemicalShift):

            # if it is not in amino acid definition
            if shift.atomLabel not in self.defined_shift_list and self.aminoAcidId != AminoAcid.find_amino_acid_id('U'):
                # warnings.warn('Attempt to add chemical shift with inconsistent label name. {} not in {} for {}.'.format(shift.atomLabel, self.defined_shift_list, AminoAcid.THREE_LETTERS_CODES[self.aminoAcidId]))
                raise Exception("Attempt to add chemical shift with inconsistent label name " + shift.atomLabel + " not in " + str(self.defined_shift_list))

            # if it is not in current list in chemical shifts
            elif shift not in self.chemicalShifts and shift.resonance is not None:
                self.chemicalShifts.append(shift)

            # if it already is in current list of chemical shifts, but resonanse doesn't agree
            elif not self.getChemicalShiftByLabel(shift.atomLabel).resonance == shift.resonance:
                raise Exception("Attempt to add chemical shift with inconsistent resonsnce: " + self.getChemicalShiftByLabel(shift.atomLabel).toString() + " vs. " + shift.toString())

    def getListOfShifts(self):
        return self.defined_shift_list

    def hasBackboneChemicalShifts(self):
        definitions = self.getListOfShifts()
        return (Atom.NITROGEN not in definitions or self.getChemicalShiftByLabel(Atom.NITROGEN) is not None) \
               and (Atom.HYDROGEN not in definitions or self.getChemicalShiftByLabel(Atom.HYDROGEN) is not None) \
               and (Atom.CARBON_ALPHA not in definitions or self.getChemicalShiftByLabel(Atom.CARBON_ALPHA) is not None) \
               and (Atom.CARBON_BETA not in definitions or self.getChemicalShiftByLabel(Atom.CARBON_BETA) is not None)

    def containsShift(self, chemicalShiftId):
        return any(shift == chemicalShiftId for shift in self.getListOfShifts())


    def currentAminoAcidDefinition(self):
        return AminoAcid.findAminoAcidDefinition(self.getAminoAcidId())

    def get_chemical_shift_completeness(self):
        complete_shifts = 0
        incomplete_shifts = 0

        for shift_label in self.defined_shift_list:
            if self.getChemicalShiftByLabel(shift_label) is not None:
                complete_shifts += 1
            else:
                incomplete_shifts += 1

        return complete_shifts, incomplete_shifts

    def isAssignmentComplete(self):
        matches = True
        for shift in self.getListOfShifts():
            # print self.AMINO_ACID_DEFINITIONS[index][i]
            # print [ elem for elem in self.chemicalShifts if elem.atomLabel == self.AMINO_ACID_DEFINITIONS[index][i] ]
            matches &= len([elem for elem in self.chemicalShifts if elem.atomLabel == shift]) == 1

        return matches

    def toString(self, full=False):
        output = self.getThreeLettersCode() + ", "

        if full:
            if len(self.chemicalShifts) == 0:
                output += "empty shift list"
            else:
                output += "Shifts (sorted): "
                output += str.join(' ', sorted([shift.atomLabel + '(' + str(shift.resonance) + ')' for shift in self.chemicalShifts]))

        return output
