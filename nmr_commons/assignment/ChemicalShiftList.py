import os
import pickle

from nmr_commons.sequence.Atom import Atom
from nmr_commons.sequence.ProteinSequence import AminoAcid
from nmr_commons.sequence.ProteinSequence import ProteinSequence
from nmr_commons.assignment.ChemicalShift import ChemicalShift


class ChemicalShiftList(ProteinSequence):
    # Supported formats
    TALOS_FORMAT = "TALOS"
    SPARKY_FORMAT = "SPARKY"
    CSV_FORMAT = "CSV"
    PICKLE_FORMAT = "PICKLE"

    def __init__(self):
        super(ChemicalShiftList, self).__init__()

    #
    #                                            SHORTCUTS
    #
    def addChemicalShift(self, chemicalShift, position):
        self.getAminoAcidById(position).addChemicalShift(chemicalShift)

    def getChemicalShifts(self, position):
        return self.getAminoAcidById(position).getChemicalShifts()

    #
    #                                           IO METHODS
    #

    def deserialize_z(self, z, sequence, listOfConsideredShifts):

        id = 0
        for aminoAcid in sequence:

            aa = AminoAcid(aminoAcid.getThreeLettersCode())
            self.addAminoAcid(aa)
            list = aminoAcid.getAllShiftNamesFromDefinition()

            for elem in list:

                if Atom.getAtomTypeFromLabel(elem) != Atom.ATOM_UNKNOWN and (elem in listOfConsideredShifts):

                    shift = z[id]
                    if not shift == None:
                        aa.addChemicalShift(ChemicalShift(elem[0], elem, shift, -1))
                    id = id + 1

    def serialize_z(self, listOfConsideredShifts):
        z = []
        for aa in self.getSequence():
            z = z + aa.serialize_z(listOfConsideredShifts)
        return z

    def getShiftOrder(self, listOfConsideredShifts):
        z = []
        for aa in self.getSequence():
            z = z + aa.getShiftOrder()
        return [elem for elem in z if elem in listOfConsideredShifts]

    def get_shifts_number(self):
        number = 0
        for aa in self.getSequence():
            number += len(aa.chemicalShifts)
        return number

    def saveAs(self, path, type=CSV_FORMAT):

        if type == ChemicalShiftList.CSV_FORMAT:
            self.saveAsCSV(path)
        elif type == ChemicalShiftList.PICKLE_FORMAT:
            pickle.dump(self.sequence, open(path, "wb"))
        elif type == ChemicalShiftList.TALOS_FORMAT:
            self.saveAsTAB(path)

    def get_chemical_shift_completeness(self):
        complete_shifts, incomplete_shifts = zip(*[aminoAcid.get_chemical_shift_completeness()
                                                   for aminoAcid in self.sequence])
        ##########

        # defi = aminoAcid.getAllShiftNamesFromAminoAcidDefinition(aminoAcid.currentAminoAcidDefinition())
        # print(aminoAcid.getThreeLettersCode())
        # print("Teoretyczne: " + str(defi))
        # print("Rzeczywiste: " + str([(x.atomLabel, x.resonance)for x in aminoAcid.getChemicalShifts()]))
        # print(result)

        ###########

        return float(sum(complete_shifts)) / sum(complete_shifts + incomplete_shifts)

    def saveAsCSV(self, path):
        f = open(path, 'w')
        f.write("AminoAcidId, AminoAcidLabel, AtomType, AtomLabel, Resonance, AmbiguityCode\n")

        for aminoAcidId in range(0, len(self.sequence)):

            aminoAcid = self.sequence[aminoAcidId]

            for shift in aminoAcid.getChemicalShifts():
                f.write(str(aminoAcidId) + "," + aminoAcid.getThreeLettersCode() + ","
                        + shift.atomType + "," + shift.atomLabel + "," + str(shift.resonance) + "," + str(
                    shift.ambiguityCode) + "\n")
        f.close()

    def saveAsTAB(self, path):
        f = open(path, 'w')

        for aminoAcidId in range(0, len(self.sequence)):
            if aminoAcidId % 50 == 0:
                f.write("\nDATA SEQUENCE")
            if aminoAcidId % 10 == 0:
                f.write(" ")

            aminoAcid = self.sequence[aminoAcidId]
            f.write(aminoAcid.getOneLetterCode())

        f.write("\n\nVARS RESID RESNAME ATOMNAME SHIFT \nFORMAT %4d %1s %4s %8.3f\n\n")

        talos_shifts = [Atom.NITROGEN, Atom.HYDROGEN_ALPHA, Atom.HYDROGEN_ALPHA2, Atom.HYDROGEN_ALPHA3, Atom.CARBON,
                        Atom.CARBON_ALPHA, Atom.CARBON_BETA]

        for i in range(0, len(self.sequence)):
            aminoAcid = self.sequence[i]
            for shift in aminoAcid.chemicalShifts:
                if shift.atomLabel in talos_shifts:
                    f.write("{:>4} {:>1} {:>4} {:>10.3f}\n".format(i + 1, aminoAcid.getOneLetterCode(), shift.atomLabel,
                                                                   shift.resonance))
                elif shift.atomLabel == Atom.HYDROGEN:
                    f.write("{:>4} {:>1} {:>4} {:>10.3f}\n".format(i + 1, aminoAcid.getOneLetterCode(), "HN",
                                                                   shift.resonance))

        f.close()

    def save_as_prot(self, path):
        with open(path, 'w') as f:
            i = 1
            for amino_acid_id in range(0, len(self.sequence)):
                occured = []
                for chemical_shift in self.sequence[amino_acid_id].chemicalShifts:
                    if chemical_shift.atomLabel not in occured:
                        f.write('{:>4} {:7.3f} {:7.3f} {:<5} {:>4}\n'.format(i, chemical_shift.resonance, 0,
                                                                         chemical_shift.atomLabel, amino_acid_id + 1))
                        occured += [chemical_shift.atomLabel]
                        i += 1

    # Methods return fragment of chemical shift list as chemical shift list
    def getFragmentOfChemmicalShiftList(self, start_Idx, end_Idx):
        newChemicalShiftList = ChemicalShiftList()

        for id in range(start_Idx, end_Idx + 1):
            newChemicalShiftList.addAminoAcid(self.getAminoAcidById(id))

        return newChemicalShiftList

    # Method finds the longest chain of amino acids within the structure that has fully assigned backbone
    # It was developed because very often at least one amino acid in not assigned in BMRB (usually the first several amino acids)
    # Mathod returns (start_Idx, end_Idx, length) where start_Idx is an index of first fully assigned amino acid and end_Idx the last one
    def getMaximumFullyAssignedFragment(self):

        start_Idx = 0
        end_Idx = 0
        length = 0

        tmp_start_Idx = 0
        tmp_end_Idx = 0
        tmp_length = 0
        for id in range(0, self.get_number_of_amino_acids()):

            # if backbone of amino acid is fully anotated
            if self.getAminoAcidById(id).hasBackboneChemicalShifts():

                # update tmp values
                tmp_start_Idx = tmp_start_Idx if tmp_length > 0 else id
                tmp_end_Idx = id
                tmp_length = tmp_end_Idx - tmp_start_Idx + 1

                if tmp_length > length:
                    length = tmp_length
                    start_Idx = tmp_start_Idx
                    end_Idx = tmp_end_Idx
            else:
                tmp_length = 0

        return [start_Idx, end_Idx, length]

    def isBackboneAssignmentComplete(self):

        for aminoAcid in self.sequence:
            if not aminoAcid.hasBackboneChemicalShifts():
                print(aminoAcid.toString(full=True))
                return False

        return True

    def loadFromFile(self, path, type=CSV_FORMAT):
        filename, file_extension = os.path.splitext(path)

        if type == ChemicalShiftList.CSV_FORMAT:
            self.loadFromCSVFile(path)
        elif type == ChemicalShiftList.PICKLE_FORMAT:
            f = open(path, "rb")
            self.sequence = pickle.load(f)
            f.close()
        elif type == ChemicalShiftList.SPARKY_FORMAT:
            self.loadFromSparkyFile(path)
        else:
            raise Exception("Unsupported chemical shift format: " + file_extension)

    def loadFromSparkyFile(self, path):
        f = open(path, 'r')

        # skip header
        f.readline()
        f.readline()

        # read file
        for line in f:
            # i.e. ['D108', 'CA', '13C', '54.342', '0.066', '4\n']
            fields = list(filter(lambda x: x != '', line.split(" ")))

            aminoAcidID = fields[0][1:]
            atomType = ''.join([i for i in fields[2] if not i.isdigit()])
            atomLabel = fields[1] if fields[1] != "HN" else "H"
            resonance = float(fields[3])
            ambiguityCode = -1

            if aminoAcidID.isdigit():
                self.addChemicalShift(ChemicalShift(atomType, atomLabel, resonance, ambiguityCode),
                                      int(aminoAcidID) - 1)

        # close file
        f.close()

    def loadFromCSVFile(self, path):
        f = open(path, 'r')

        # skip header
        f.readline()

        # read file
        for line in f:
            # AminoAcidId, AminoAcidLabel, AtomType, AtomLabel, w1, w2, w3
            fields = line.split(",")
            print(fields)
            self.addChemicalShift(ChemicalShift(fields[2], fields[3], float(fields[4]), float(fields[5])),
                                  int(fields[0]))

        # close file
        f.close()

    def load_from_prot(self, path):
        with open(path, 'r') as f:
            for line in f:
                fields = line.split()

                amino_id = int(fields[4])-1
                atom_label = fields[3]
                if atom_label.startswith('Q'):
                    atom_labels = Atom.get_atom_name_for_pseudo_atom(atom_label, self.getAminoAcidById(amino_id))
                else:
                    atom_labels = [atom_label]
                # atom_type = Atom.getAtomTypeFromLabel(atom_label)
                resonance = float(fields[1])
                ambiguity_code = -1

                if resonance != 999.0:
                    for label in atom_labels:
                        atom_type = Atom.getAtomTypeFromLabel(label)
                        self.addChemicalShift(ChemicalShift(atom_type, label, resonance, ambiguity_code), amino_id)

    def toString(self, sequenceOnly=False):

        if sequenceOnly:
            return super(ChemicalShiftList, self).toString()
        else:
            result = "=== List of chemical shifts ===\n"
            idNum = 1
            for aminoAcid in self.sequence:
                result += str(idNum) + ":\n"
                result += aminoAcid.toString(True) + "\n"
                idNum += 1
            return result
