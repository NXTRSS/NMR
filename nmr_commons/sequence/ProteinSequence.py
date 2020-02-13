import os

from nmr_commons.sequence.AminoAcid import AminoAcid


class ProteinSequence(object):

    def __init__(self):
        self.sequence = []

    def addAminoacidsFromOneLetterCodes(self, codeStr):

        for letter in codeStr:
            self.addAminoAcid(AminoAcid(letter))

    def get_weight_in_Da(self):
        return reduce( lambda acc, x: acc + x.get_weight_in_Da(), self.sequence, 0)

    # Return number of unknown amino acids. Unknown amino acids may appear in a sequence for instance when *.str file
    # with improper one- threeletter code is loaded
    def get_number_of_unknown_amino_acids(self):

        unknown_elems = 0

        for elem in self.sequence:
            if elem.definition == AminoAcid.UNKNOWN:
                unknown_elems += 1

        return unknown_elems

    def addAminoAcid(self, aminoAcid):
        if isinstance(aminoAcid, AminoAcid):
            self.sequence.append(aminoAcid)
        else:
            raise Exception("Attempt to insert invalid amino acid object")

    def get_number_of_amino_acids(self):
        return len(self.sequence)

    def getSequence(self):
        return self.sequence

    def getAminoAcidById(self, id):
        return self.sequence[id] if 0 <= id < len(self.sequence) else None

    #                                           IO METHODS
    def saveSequence(self, path):
        f = open(path, 'w')
        f.write("> Protein sequence\n")
        for elem in range(0, len(self.sequence)):
            f.write(self.getAminoAcidById(elem).getOneLetterCode())
        f.close()

    def save_sequence_with_three_letters_codes(self, path):
        with open(path, 'w') as f:
            for elem in range(0, len(self.sequence)):
                f.write(self.getAminoAcidById(elem).getThreeLettersCode() + '\n')

    def loadSequence(self, path):
        filename, fileExtension = os.path.splitext(path)
        if fileExtension == '.fasta':
            self.load_from_fasta_file(path)
        else:
            self.load_sequence_from_file(path)
            # raise "Unsupported format while loading protein sequence (supported formats *.fasta)"

    def load_sequence_from_file(self, path):
        with open(path, 'r') as f:
            for line in f:
                identifier = ''.join([i for i in line if not i.isdigit()])
                identifier = identifier.strip()
                self.addAminoAcid(AminoAcid(identifier))

    def load_from_fasta_file(self, path):
        f = open(path, 'r')
        f.readline()
        lines = f.readlines()
        for line in lines:
            for ch in line:
                if not ch.isspace():
                    self.addAminoAcid(AminoAcid(ch))
        f.close()

    def toString(self):
        result = ""
        for elem in range(0, len(self.sequence)):
            result += str(elem) + ":" + self.sequence[elem].toString() + " "
        return "Pritein sequence object is empty" if len(result) == 0 else result

