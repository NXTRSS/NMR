class Atom:
    # ATOM TYPES
    ATOM_HYDROGEN = "H"
    ATOM_CARBON = "C"
    ATOM_CO = "CO"
    ATOM_NITROGEN = "N"
    ATOM_UNKNOWN = "U"

    # ATOM LABELS
    HYDROGEN = "H"
    HYDROGEN_ALPHA = "HA"
    HYDROGEN_ALPHA2 = "HA2"
    HYDROGEN_ALPHA3 = "HA3"
    HYDROGEN_BETA = "HB"
    HYDROGEN_BETA2 = "HB2"
    HYDROGEN_BETA3 = "HB3"
    HYDROGEN_GAMMA = "HG"
    HYDROGEN_GAMMA1 = "HG1"
    HYDROGEN_GAMMA12 = "HG12"
    HYDROGEN_GAMMA13 = "HG13"
    HYDROGEN_GAMMA2 = "HG2"
    HYDROGEN_GAMMA3 = "HG3"
    HYDROGEN_DELTA1 = "HD1"
    HYDROGEN_DELTA2 = "HD2"
    HYDROGEN_DELTA3 = "HD3"
    HYDROGEN_EPSILON = "HE"
    HYDROGEN_EPSILON1 = "HE1"
    HYDROGEN_EPSILON2 = "HE2"
    HYDROGEN_EPSILON3 = "HE3"
    HYDROGEN_HH = "HH"
    HYDROGEN_H11 = "HH11"
    HYDROGEN_H12 = "HH12"
    HYDROGEN_HH2 = "HH2"
    HYDROGEN_H21 = "HH21"
    HYDROGEN_H22 = "HH22"
    HYDROGEN_DELTA21 = "HD21"
    HYDROGEN_DELTA22 = "HD22"
    HYDROGEN_EPSILON21 = "HE21"
    HYDROGEN_EPSILON22 = "HE22"
    HYDROGEN_ZETA = "HZ"
    HYDROGEN_ZETA2 = "HZ2"
    HYDROGEN_ZETA3 = "HZ3"
    ALL_HYDROGENS = [HYDROGEN, HYDROGEN_ALPHA, HYDROGEN_ALPHA2, HYDROGEN_ALPHA3, HYDROGEN_BETA, HYDROGEN_BETA2,
                     HYDROGEN_BETA3,
                     HYDROGEN_GAMMA, HYDROGEN_GAMMA1, HYDROGEN_GAMMA12, HYDROGEN_GAMMA13, HYDROGEN_GAMMA2,
                     HYDROGEN_GAMMA3,
                     HYDROGEN_DELTA1, HYDROGEN_DELTA2, HYDROGEN_DELTA3, HYDROGEN_EPSILON, HYDROGEN_EPSILON1,
                     HYDROGEN_EPSILON2,
                     HYDROGEN_EPSILON3, HYDROGEN_HH, HYDROGEN_H11, HYDROGEN_H12, HYDROGEN_HH2, HYDROGEN_H21,
                     HYDROGEN_H22,
                     HYDROGEN_DELTA21, HYDROGEN_DELTA22, HYDROGEN_EPSILON21, HYDROGEN_EPSILON22, HYDROGEN_ZETA,
                     HYDROGEN_ZETA2, HYDROGEN_ZETA3]

    CARBON = "C"
    CARBON_ALPHA = "CA"
    CARBON_BETA = "CB"
    CARBON_GAMMA = "CG"
    CARBON_GAMMA1 = "CG1"
    CARBON_GAMMA2 = "CG2"
    CARBON_DELTA = "CD"
    CARBON_DELTA1 = "CD1"
    CARBON_DELTA2 = "CD2"
    CARBON_EPSILON = "CE"
    CARBON_EPSILON1 = "CE1"
    CARBON_EPSILON2 = "CE2"
    CARBON_EPSILON3 = "CE3"
    CARBON_ZETA = "CZ"
    CARBON_ZETA1 = "CZ1"
    CARBON_H2 = "CH2"
    CARBON_ZETA2 = "CZ2"
    CARBON_ZETA3 = "CZ3"
    ALL_CARBONS = [CARBON, CARBON_ALPHA, CARBON_BETA, CARBON_GAMMA, CARBON_GAMMA1, CARBON_GAMMA2, CARBON_DELTA,
                   CARBON_DELTA1, CARBON_DELTA2, CARBON_EPSILON,
                   CARBON_EPSILON1, CARBON_EPSILON2, CARBON_EPSILON3, CARBON_ZETA, CARBON_ZETA1, CARBON_H2,
                   CARBON_ZETA2, CARBON_ZETA3]

    NITROGEN = "N"
    NITROGEN_EPSILON = "NE"
    NITROGEN_EPSILON2 = "NE2"
    NITROGEN_NH1 = "NH1"
    NITROGEN_NH2 = "NH2"
    NITROGEN_DELTA2 = "ND2"
    NITROGEN_ZETA = "NZ"
    NITROGEN_EPSILON1 = "NE1"
    NITROGEN_DELTA1 = "ND1"
    ALL_NITROGENS = [NITROGEN, NITROGEN_EPSILON, NITROGEN_EPSILON2, NITROGEN_NH1, NITROGEN_NH2, NITROGEN_DELTA2,
                     NITROGEN_ZETA,
                     NITROGEN_EPSILON1, NITROGEN_DELTA1]

    OXYGEN_DELTA = "OD"
    OXYGEN_GAMMA = "OG"
    OXYGEN_GAMMA1 = "OG1"
    OXYGEN_EPSILON = "OE"
    OXYGEN_ETA = "OH"
    SULFUR_EPSILON = "SE"
    SULFUR_GAMMA = "SG"
    SULFUR_DELTA = "SD"
    ALL_SHIFTS = ALL_HYDROGENS + ALL_CARBONS + ALL_NITROGENS

    @staticmethod
    def get_atom_weight(label):
        if label == 'NA':
            return 100
        elif label in ['N', 'NH', '15N']:
            return 3
        elif label in ['C', 'CH', '13C', 'c']:
            return 2
        elif label in ['H', 'HN', 'HC', '1H', 'h']:
            return 1
        else:
            return -1

    @staticmethod
    def get_universal_atom_label(label):
        if label in ['N', 'NH', '15N', 'N1', 'N2']:
            return 'N'
        elif label in ['C', 'CH', '13C', 'C1', 'C2', 'C3', 'CO', 'CA', 'CB']:
            return 'C'
        elif label in ['H', 'HN', 'HC', '1H', 'H1', 'H2']:
            return 'H'
        else:
            return label

    @staticmethod
    def getAtomTypeFromLabel(atomLabel):
        if "H" in atomLabel and atomLabel.index("H") == 0:
            return Atom.ATOM_HYDROGEN
        if "C" in atomLabel and len(atomLabel) == 1:
            return Atom.ATOM_CO
        if "C" in atomLabel and atomLabel.index("C") == 0:
            return Atom.ATOM_CARBON
        if "N" in atomLabel and atomLabel.index("N") == 0:
            return Atom.ATOM_NITROGEN
        return Atom.ATOM_UNKNOWN

    @staticmethod
    def get_atom_name_for_pseudo_atom(pseudo_atom_label, amino_acid):
        if pseudo_atom_label == 'QQD':
            return Atom.get_atom_name_for_pseudo_atom('QD1', amino_acid) + \
                   Atom.get_atom_name_for_pseudo_atom('QD2', amino_acid)
        elif pseudo_atom_label == 'QQG':
            return Atom.get_atom_name_for_pseudo_atom('QG1', amino_acid) + \
                   Atom.get_atom_name_for_pseudo_atom('QG2', amino_acid)
        elif pseudo_atom_label == 'QR':
            return Atom.get_atom_name_for_pseudo_atom('QD', amino_acid) + \
                   Atom.get_atom_name_for_pseudo_atom('QE', amino_acid) + \
                   Atom.get_atom_name_for_pseudo_atom('HZ', amino_acid)

        else:
            replaced_label = pseudo_atom_label.replace('Q', 'H')
            matching_labels = [label for label in amino_acid.defined_shift_list if label.startswith(replaced_label)]
            result = []
            for label in matching_labels:
                if not amino_acid.getChemicalShiftByLabel(label):
                    result += [label]
            return result
