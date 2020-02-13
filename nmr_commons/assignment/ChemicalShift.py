class ChemicalShift:

    def __init__(self, atomType, atomLabel, resonance, ambiguityCode):
        self.atomType = atomType  # N, C, H
        self.atomLabel = atomLabel  # CA, CB, CG
        self.resonance = float(resonance)
        self.ambiguityCode = ambiguityCode

    def toString(self):
        return self.atomType + " " + self.atomLabel + " " + str(self.resonance) + " " + str(self.ambiguityCode)



