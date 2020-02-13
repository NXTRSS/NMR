import pickle
import re
from timeit import default_timer as timer
from nmr_commons.assignment.ChemicalShift import ChemicalShift
from nmr_commons.assignment.ChemicalShiftList import ChemicalShiftList
from nmr_commons.sequence.ProteinSequence import *
import nmr_commons.database.bmrb as bmrb
from nmr_commons.utils.ConsoleUtils import *
from nmr_environment import Settings


class BMRBDatabase:

    def __init__(self, bmrb_root, max_num_files=20):
        self.parsed_db = []
        self.bmrb_root = ''

        if os.path.isdir(bmrb_root):
            self.bmrb_root = bmrb_root
        else:
            start = timer()
            self.parsed_db = pickle.load(open(bmrb_root, 'rb'))
            print('Reading database of {} took {:.2f}s'.format(len(self.parsed_db), timer() - start))

        self.max_num_files = max_num_files

    def save(self, path):
        pickle.dump(self.parsed_db, open(path, 'wb'))

    @staticmethod
    def get_all_bmrb_ids():
        return [int(re.findall("\d+", f)[0]) for f in os.listdir(Settings.BMRB_ROOT)]

    def parse_bmrb_db(self):
        self.parsed_db = []
        bmrb_ids = self.get_train_ids()
        num_processed = 0
        num_invalid = 0
        start = timer()
        num_to_process = min(self.max_num_files, len(bmrb_ids))

        print('Number of files to process: {}'.format(num_to_process))

        for bmrb_id in bmrb_ids[:num_to_process]:
            try:
                self.parsed_db.append(self.parse_bmrb_file('{}bmr{}.str'.format(self.bmrb_root, bmrb_id)))
                num_processed += 1
            except:
                num_invalid += 1
                print('WARNING: File {} cannot be processed because of invalid content.'.format(bmrb_id))

            # report progress
            printProgress(num_processed, num_to_process)

        print('\nTime elapsed: {:.2f}s'.format(timer() - start))
        print('Invalid files: {}, valid files: {}'.format(num_invalid, num_processed))

    @staticmethod
    def get_train_ids():
        with open("../../chemical_shift_assignment/data/TrainIDx.txt") as f:
            ids = [int(line) for line in f]
        return ids

    @staticmethod
    def getTestIdx():
        f = open("../Database/TestIDx.txt")
        listTr = [int(line) for line in f]
        f.close()
        return listTr

    @staticmethod
    def get_related_pdb_entries(filePath):
        bmrbEntry = bmrb.entry.fromFile(filePath)
        relatedEntries = bmrbEntry.getLoopsByCategory("Related_entries")

        for entry in relatedEntries:

            if entry['Database_name'][0].lower() == 'pdb':
                return entry['Database_accession_code'][0]

        return None

    #
    #                                           Parse single bmrb file
    #
    @staticmethod
    def parse_bmrb_file(filePath):
        bmrbEntry = bmrb.entry.fromFile(filePath)

        # Chemical shift list object
        chemicalShifts = ChemicalShiftList()

        # Read protein sequence
        seq = bmrbEntry.getSaveframesByCategory("entity")
        seq = str(seq[0]["Polymer_seq_one_letter_code"][0].replace('\n', ''))

        if len(seq) == 0:
            raise Exception("Can not retrieve sequence information from file: " + filePath)

        chemicalShifts.addAminoacidsFromOneLetterCodes(seq)

        # Read chemical shifts
        bmrbEntryData = bmrbEntry.getSaveframesByCategory('assigned_chemical_shifts')[0][
            '_Atom_chem_shift'].getDataByTag(['Seq_ID', 'Comp_ID', 'Atom_type', 'Atom_ID', 'Val', 'Ambiguity_code'])

        for entry in bmrbEntryData:
            index = int(entry[0]) - 1
            aminoAcid = chemicalShifts.getSequence()[index]
            chemShift = ChemicalShift(*entry[2:6])

            # in case shift H*1=Y is saved in BMRB as H*11=Y, H*12=Y, H*13=Y. These two options are possible if three hyydrogens give th same chemical shift
            if chemShift.atomType == 'H' and len(chemShift.atomLabel) == 4 and chemShift.atomLabel not in aminoAcid.defined_shift_list:
                chemShift.atomLabel = chemShift.atomLabel[:-1]

            # in case shift H* is saved as H*1 or H*2
            if chemShift.atomType == 'H' and len(chemShift.atomLabel) == 3 and chemShift.atomLabel not in aminoAcid.defined_shift_list:
                chemShift.atomLabel = chemShift.atomLabel[:-1]

            if AminoAcid.find_amino_acid_id(entry[1]) == chemicalShifts.getAminoAcidById(index).aminoAcidId:
                chemicalShifts.addChemicalShift(chemShift, int(entry[0])-1)
            else:
                return None

        return chemicalShifts