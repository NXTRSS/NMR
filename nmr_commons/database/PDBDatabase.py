from nmr_commons.structure.SecondaryStructure import SecondaryStructure
from nmr_environment import Settings


class PDBDatabase:

    @staticmethod
    def get_secondary_structure(entry_pdb_code):
        key = ">" + entry_pdb_code

        f = open(Settings.PDB_SECONDARY_STRUCTURES)
        line = f.readline()

        secondary_structure = ""

        while len(line) > 0:

            # Check it it is proper secondary structure section
            if line.startswith(key):
                if line.endswith('secstr\n'):

                    # Read secondary structure
                    while True:
                        structure_line = f.readline()
                        if not structure_line.startswith(">"):
                            secondary_structure += structure_line
                        else:
                            break
                    break

            line = f.readline()
        f.close()

        secondary_structure = secondary_structure.replace('\n','');
        return SecondaryStructure(secondary_structure) if len(secondary_structure) > 0 else None