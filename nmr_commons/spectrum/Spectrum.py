from __future__ import print_function

import warnings

import nmrglue as ng
import numpy as np
# from scipy.signal import fftconvolve
from nmr_commons.sequence.AminoAcid import AminoAcid
from nmr_commons.sequence.Atom import Atom
import scipy.ndimage as ndimage
from nmr_commons.peak_picking.PeakList import PeakList
from nmr_commons.spectrum.Data import Data
from nmr_environment import Settings
from nmr_commons.spectrum.Bundle import Bundle
import matplotlib.pyplot as plt
import skimage.measure as mes

from multiprocessing import Pool
from timeit import default_timer as timer
from skimage import exposure

_PARALLEL_CONTEXT = None


def _parallel_work_impl(elem):

    obj, dx, dy, dz, dxo, dyo, dzo, normalization_mode = _PARALLEL_CONTEXT
    if normalization_mode == Spectrum.NORMALIZATION_EQ_0_255:
        tmp = obj.get_cube(elem[0], elem[1], elem[2], dx, dy, dz, obj.eq_hist_data)
    else:
        tmp = obj.get_cube(elem[0], elem[1], elem[2], dx, dy, dz)

    # Minimap 1
    # tmp2 = obj.get_cube(elem[0], elem[1], int(obj.get_size()[2]/2.0), dx, dy, int(obj.get_size()[2]/2.0))
    # tmp2 = (np.abs(tmp2) - obj.median_abs) / obj.std_abs
    # tmp2 = mes.block_reduce(tmp2, block_size=(1, 1, int(obj.get_size()[2]/dzo)), func=np.max)
    # tmp2 = ndimage.zoom(tmp2, (1.0 * dxo/tmp2.shape[0], 1.0 * dyo/tmp2.shape[1], 1.0 * dzo/tmp2.shape[2]), order=1)


    # Minimap 2
    # tmp3 = obj.get_cube(elem[0], int(obj.get_size()[1]/2.0), elem[2], dx, int(obj.get_size()[1]/2.0), dz)
    # tmp3 = (np.abs(tmp3) - obj.median_abs) / obj.std_abs
    # tmp3 = mes.block_reduce(tmp3, block_size=(1, int(obj.get_size()[1]/dzo), 1), func=np.max)
    # tmp3 = ndimage.zoom(tmp3, (1.0 * dxo/tmp3.shape[0], 1.0 * dyo/tmp3.shape[1], 1.0 * dzo/tmp3.shape[2]), order=1)

    # plt.figure()
    # print(elem)
    # layer = (np.abs(obj.data[elem[0], :, :]) - obj.median_abs) / obj.std_abs
    # plt.imshow(np.squeeze(layer))
    # plt.clim(0,0.3)
    # plt.figure()
    # plt.imshow(np.squeeze(tmp2))
    # plt.clim(0,3)
    # plt.figure()
    # plt.imshow(np.squeeze(tmp3))
    # plt.clim(0,3)
    # plt.show()


    if tmp.shape != (dxo, dyo, dzo):
        tmp = ndimage.zoom(tmp, (1.0 * dxo/(2 * dx + 1), 1.0 * dyo/(2 * dy + 1), 1.0 * dzo/(2 * dz + 1)), order=1)

    if normalization_mode == Spectrum.NORMALIZATION_GAUSS:
        tmp = (np.abs(tmp) - obj.median_abs) / obj.std_abs
    elif normalization_mode == Spectrum.NORMALIZATION_LAPLACE:
        tmp = (np.abs(tmp) - obj.median_abs) / obj.mad
    elif normalization_mode == Spectrum.NORMALIZATION_0_255:
        tmp = np.abs(tmp)
        tmp[tmp < 1] = 1
        tmp = np.log(tmp)
        tmp = tmp / obj.log_mean
        tmp = tmp / (np.log(obj.max_absolute_value) / obj.log_mean) * 255.0
        tmp = np.array(tmp, dtype=np.float16)
        # tmp = tmp / obj.max_absolute_value * 255.
        # tmp = (np.abs(tmp) - obj.median_abs) / obj.std_abs
        # tmp = np.abs(tmp) / ((obj.max_absolute_value - obj.median_abs) / obj.std_abs) * 255.0
    elif normalization_mode == Spectrum.NORMALIZATION_EQ_0_255:
        tmp *= 255
    elif normalization_mode == Spectrum.NORMALIZATION_HAMPEL:
        tmp = 0.5*(np.tanh(0.1 * ((tmp - obj.mean) / obj.std)) + 1)
    elif normalization_mode == Spectrum.NORMALIZATION_UNIFORM:
        tmp = np.abs(tmp) / obj.max_absolute_value
    elif normalization_mode == Spectrum.NORMALIZATION_RAYLEIGH:
        tmp = np.abs(tmp)
        tmp[tmp < 1] = 1
        tmp = np.log(tmp)
        tmp = 0.5 * np.pi * tmp / obj.log_mean
    elif normalization_mode == Spectrum.NORMALIZATION_CENTRAL_TO_ONE:
        tmp = (np.abs(tmp) - obj.median_abs) / obj.std_abs
        tmp = tmp / tmp[int(dxo/2.0), int(dyo/2.0), int(dzo/2.0)]
    elif normalization_mode == Spectrum.NORMALIZATION_NONE:
        pass
    else:
        raise Exception("Invalid normalization mode")

    # patch = np.concatenate( (tmp.transpose([1, 2, 0]), tmp2.transpose([1, 2, 0]), tmp3.transpose([1, 2, 0])), axis=2)
    patch = tmp.transpose([1, 2, 0])

    del tmp
    return patch


def _extrema_worker(layer_id):
    return _PARALLEL_CONTEXT.get_extrema_on_layer(layer_id)


# CLASS FUNCTIONS AND PROPERTIES:
#   - Defines all types of spectra
#   - Generates response maps
#   - Exports response maps to Sparky format
class Spectrum(Data):

    #
    # Axes definitions
    NITROGEN_AXIS = {'abbr': 'N', 'label': '15N', 'atom': Atom.NITROGEN, 'PPMRange': [90, 140], 'resolution': 64, 'frequency': 60}
    NITROGEN_AXIS_HD = {'label': '15N', 'atom': Atom.NITROGEN, 'PPMRange': [90, 140], 'resolution': 128, 'frequency': 60}
    CARBON_AXIS = {'abbr': 'C', 'label': '13C', 'atom': Atom.CARBON, 'PPMRange': [10, 90], 'resolution': 128, 'frequency': 150}
    CARBON_AXIS_HD = {'label': '13C', 'atom': Atom.CARBON, 'PPMRange': [90, 140], 'resolution': 1024, 'frequency': 60}
    CARBON_CO_AXIS = {'abbr': 'CO', 'label': '13C', 'atom': Atom.CARBON, 'PPMRange': [160, 190], 'resolution': 128, 'frequency': 150}
    HYDROGEN_AXIS = {'abbr': 'H', 'label': '1H', 'atom': Atom.HYDROGEN, 'PPMRange': [4, 12], 'resolution': 256, 'frequency': 600}
    HYDROGEN_AXIS_HD = {'label': '1H', 'atom': Atom.HYDROGEN, 'PPMRange': [4, 12], 'resolution': 256, 'frequency': 600}
    DUMMY_AXIS = {'label': 'NA', 'atom': Atom.ATOM_UNKNOWN, 'PPMRange': [0, 0], 'resolution': 1, 'frequency': -1}
    ALL_AXES = [NITROGEN_AXIS, NITROGEN_AXIS_HD, CARBON_CO_AXIS, CARBON_AXIS, CARBON_AXIS_HD, HYDROGEN_AXIS, HYDROGEN_AXIS_HD]

    #
    #
    NORMALIZATION_NONE = 0
    NORMALIZATION_GAUSS = 1
    NORMALIZATION_0_255 = 2
    NORMALIZATION_LAPLACE = 3
    NORMALIZATION_HAMPEL = 4
    NORMALIZATION_UNIFORM = 5
    NORMALIZATION_RAYLEIGH = 6
    NORMALIZATION_EQ_0_255 = 7
    NORMALIZATION_CENTRAL_TO_ONE = 8

    #
    #  Spectra definitions
    TYPE_BACKBONE = 0
    TYPE_SIDE_CHAIN = 1
    TYPE_NOESY = 2

    TYPE_THROUGH_BOND = 0
    TYPE_THROUGH_SPACE = 1

    #
    #  Spectrum parent type
    BACKBONE = 'backbone'
    SIDE_CHAIN = 'side_chain'
    NOESY = 'noesy'


    HNCA = {
        'name': "HNCA",
        'residueType': TYPE_BACKBONE,
        'experimentClass': TYPE_THROUGH_BOND,
        'peakDefinition':
            [[(0, Atom.NITROGEN), (-1, Atom.CARBON_ALPHA), (0, Atom.HYDROGEN)],
             [(0, Atom.NITROGEN), (0, Atom.CARBON_ALPHA), (0, Atom.HYDROGEN)]],
        'axes': [NITROGEN_AXIS, CARBON_AXIS, HYDROGEN_AXIS],
        'nDim': 3
    }

    HNCO = {
        'name': "HNCO",
        'residueType': TYPE_BACKBONE,
        'experimentClass': TYPE_THROUGH_BOND,
        'peakDefinition':
            [[(0, Atom.NITROGEN), (-1, Atom.CARBON), (0, Atom.HYDROGEN)]],
        'axes': [NITROGEN_AXIS, CARBON_CO_AXIS, HYDROGEN_AXIS],
        'nDim': 3
    }

    HNCACB = {
        'name': "HNCACB",
        'residueType': TYPE_BACKBONE,
        'experimentClass': TYPE_THROUGH_BOND,
        'peakDefinition':
            [[(0, Atom.NITROGEN), (-1, Atom.CARBON_ALPHA), (0, Atom.HYDROGEN)],
             [(0, Atom.NITROGEN), (-1, Atom.CARBON_BETA), (0, Atom.HYDROGEN)],
             [(0, Atom.NITROGEN), (0, Atom.CARBON_ALPHA), (0, Atom.HYDROGEN)],
             [(0, Atom.NITROGEN), (0, Atom.CARBON_BETA), (0, Atom.HYDROGEN)]],
        'axes': [NITROGEN_AXIS, CARBON_AXIS, HYDROGEN_AXIS],
        'nDim': 3
    }

    CBCACONH = {
        'name': "CBCACONH",
        'residueType': TYPE_BACKBONE,
        'experimentClass': TYPE_THROUGH_BOND,
        'peakDefinition':
            [[(0, Atom.NITROGEN), (-1, Atom.CARBON_ALPHA), (0, Atom.HYDROGEN)],
             [(0, Atom.NITROGEN), (-1, Atom.CARBON_BETA), (0, Atom.HYDROGEN)]],
        'axes': [NITROGEN_AXIS, CARBON_AXIS, HYDROGEN_AXIS],
        'nDim': 3
    }
    HNCOCBCA = CBCACONH

    HNCOCA = {
        'name': "HNCOCA",
        'residueType': TYPE_BACKBONE,
        'experimentClass': TYPE_THROUGH_BOND,
        'peakDefinition':
            [[(0, Atom.NITROGEN), (-1, Atom.CARBON_ALPHA), (0, Atom.HYDROGEN)]],
        'axes': [NITROGEN_AXIS, CARBON_AXIS, HYDROGEN_AXIS],
        'nDim': 3
    }
    CACONH = HNCOCA

    HNCACO = {
        'name': "HNCACO",
        'residueType': TYPE_BACKBONE,
        'experimentClass': TYPE_THROUGH_BOND,
        'peakDefinition':
            [[(0, Atom.NITROGEN), (-1, Atom.CARBON), (0, Atom.HYDROGEN)],
             [(0, Atom.NITROGEN), (0, Atom.CARBON), (0, Atom.HYDROGEN)]],
        'axes': [NITROGEN_AXIS, CARBON_CO_AXIS, HYDROGEN_AXIS],
        'nDim': 3
    }

    HBHACONH = {
        'name': "HBHACONH",
        'residueType': TYPE_BACKBONE,
        'experimentClass': TYPE_THROUGH_BOND,
        'peakDefinition':
            [[(0, Atom.NITROGEN), (0, Atom.HYDROGEN), (-1, Atom.HYDROGEN_ALPHA)],
             [(0, Atom.NITROGEN), (0, Atom.HYDROGEN), (-1, Atom.HYDROGEN_ALPHA2)],
             [(0, Atom.NITROGEN), (0, Atom.HYDROGEN), (-1, Atom.HYDROGEN_ALPHA3)],
             [(0, Atom.NITROGEN), (0, Atom.HYDROGEN), (-1, Atom.HYDROGEN_BETA)],
             [(0, Atom.NITROGEN), (0, Atom.HYDROGEN), (-1, Atom.HYDROGEN_BETA2)],
             [(0, Atom.NITROGEN), (0, Atom.HYDROGEN), (-1, Atom.HYDROGEN_BETA3)]],
        'axes': [NITROGEN_AXIS, HYDROGEN_AXIS, HYDROGEN_AXIS],
        'nDim': 3
    }

    NHSQC = {
        'name': "NHSQC",
        'residueType': TYPE_BACKBONE,
        'experimentClass': TYPE_THROUGH_BOND,
        'peakDefinition':
            [[(0, Atom.NITROGEN), (0, Atom.HYDROGEN)]],
        'axes': [NITROGEN_AXIS_HD, HYDROGEN_AXIS_HD],
        'nDim': 2
    }

    CHSQC = {
        'name': "CHSQC",
        'residueType': TYPE_SIDE_CHAIN,
        'experimentClass': TYPE_THROUGH_BOND,
        'mainAtomTypes': [Atom.ATOM_CARBON],
        'hydroBondAtomTypes': [Atom.ATOM_CARBON],
        'neighborhood': 0,
        'axes': [CARBON_AXIS_HD, HYDROGEN_AXIS_HD],
        'nDim': 2
    }

    HCCHTOCSY = {
        'name': "HCCHTOCSY",
        'residueType': TYPE_SIDE_CHAIN,
        'experimentClass': TYPE_THROUGH_BOND,
        'mainAtomTypes': [Atom.ATOM_CARBON],
        'hydroBondAtomTypes': [Atom.ATOM_CARBON],
        'neighborhood': np.inf,
        'axes': [CARBON_AXIS, HYDROGEN_AXIS, HYDROGEN_AXIS],
        'nDim': 3
    }

    HCCHCOSY = {
        'name': "HCCHCOSY",
        'residueType': TYPE_SIDE_CHAIN,
        'experimentClass': TYPE_THROUGH_BOND,
        'mainAtomTypes': [Atom.ATOM_CARBON],
        'hydroBondAtomTypes': [Atom.ATOM_CARBON],
        'neighborhood': 1,
        'axes': [CARBON_AXIS, HYDROGEN_AXIS, HYDROGEN_AXIS],
        'nDim': 3
    }

    CCCONH = {
        'name': "CCCONH",
        'residueType': TYPE_SIDE_CHAIN,
        'experimentClass': TYPE_THROUGH_BOND,
        'mainAtomTypes': [Atom.ATOM_NITROGEN],
        'hydroBondAtomTypes': [Atom.ATOM_CARBON],
        'neighborhood': np.inf,
        'axes': [NITROGEN_AXIS, HYDROGEN_AXIS, CARBON_AXIS],
        'shifts_in_amino_acid_chain': [-1],
        'nDim': 3
    }

    HCCONH = {
        'name': "HCCONH",
        'residueType': TYPE_SIDE_CHAIN,
        'experimentClass': TYPE_THROUGH_BOND,
        'mainAtomTypes': [Atom.ATOM_NITROGEN],
        'hydroBondAtomTypes': [Atom.ATOM_CARBON],
        'neighborhood': np.inf,
        'axes': [NITROGEN_AXIS, HYDROGEN_AXIS, HYDROGEN_AXIS],
        'shifts_in_amino_acid_chain': [-1],
        'nDim': 3
    }

    N15TOCSY = {
        'name': "N15TOCSY",
        'residueType': TYPE_SIDE_CHAIN,
        'experimentClass': TYPE_THROUGH_BOND,
        'mainAtomTypes': [Atom.ATOM_NITROGEN],
        'hydroBondAtomTypes': [Atom.ATOM_NITROGEN, Atom.ATOM_CARBON],
        'neighborhood': np.inf,
        'axes': [NITROGEN_AXIS, HYDROGEN_AXIS, HYDROGEN_AXIS],
        'nDim': 3
    }

    HHNOESY = {
        'name': "HHNOESY",
        'residueType': TYPE_SIDE_CHAIN,
        'experimentClass': TYPE_THROUGH_SPACE,
        'mainAtomTypes': [Atom.ATOM_NITROGEN, Atom.ATOM_CARBON],
        'hydroBondAtomTypes': [Atom.ATOM_NITROGEN, Atom.ATOM_CARBON],
        'neighborhood': np.inf,
        'shifts_in_amino_acid_chain': [-1, 0],
        'axes': [HYDROGEN_AXIS, HYDROGEN_AXIS],
        'nDim': 2
    }

    N15NOESY = {
        'name': "N15NOESY",
        'residueType': TYPE_SIDE_CHAIN,
        'experimentClass': TYPE_THROUGH_SPACE,
        'mainAtomTypes': [Atom.ATOM_NITROGEN],
        'hydroBondAtomTypes': [Atom.ATOM_NITROGEN, Atom.ATOM_CARBON],
        'neighborhood': np.inf,
        'shifts_in_amino_acid_chain': [-1, 0],
        'axes': [NITROGEN_AXIS, HYDROGEN_AXIS, HYDROGEN_AXIS],
        'nDim': 3
    }

    C13NOESY = {
        'name': "C13NOESY",
        'residueType': TYPE_SIDE_CHAIN,
        'experimentClass': TYPE_THROUGH_SPACE,
        'mainAtomTypes': [Atom.ATOM_CARBON],
        'hydroBondAtomTypes': [Atom.ATOM_NITROGEN, Atom.ATOM_CARBON],
        'neighborhood': np.inf,
        'shifts_in_amino_acid_chain': [-1, 0],
        'axes': [CARBON_AXIS, HYDROGEN_AXIS, HYDROGEN_AXIS],
        'nDim': 3
    }

    NAME_DEF_DICT = {
        HNCA['name']:         HNCA,
        HNCO['name']:         HNCO,
        HNCACB['name']:       HNCACB,
        CBCACONH['name']:     CBCACONH,
        CACONH['name']:       CACONH,
        HNCACO['name']:       HNCACO,
        NHSQC['name']:        NHSQC,
        CHSQC['name']:        CHSQC,
        HBHACONH['name']:     HBHACONH,
        HCCHTOCSY['name']:    HCCHTOCSY,
        HCCHCOSY['name']:     HCCHCOSY,
        N15TOCSY['name']:     N15TOCSY,
        HHNOESY['name']:      HHNOESY,
        N15NOESY['name']:     N15NOESY,
        C13NOESY['name']:     C13NOESY,
        HCCONH['name']:       HCCONH,
        CCCONH['name']:       CCCONH
    }

    CYANA_TYPES_DICT = {
        HNCA['name']:        'HNCA',
        HNCO['name']:        'HNCO',
        HNCACB['name']:      'CBCANH',
        CBCACONH['name']:    'CBCAcoNH',
        CACONH['name']:      'HNcoCA',
        HNCACO['name']:      'HNcaCO',
        NHSQC['name']:       'N15HSQC',
        CHSQC['name']:       'C13HSQC',
        HBHACONH['name']:    'HBHAcoNH',
        HCCHTOCSY['name']:   'HCCHTOCSY',
        HCCHCOSY['name']:    'HCCHCOSY',
        N15TOCSY['name']:    'N15TOCSY',
        HHNOESY['name']:     'NOESY',
        N15NOESY['name']:    'N15NOESY',
        C13NOESY['name']:    'C13NOESY',
        HCCONH['name']:      'HCcoNH',
        CCCONH['name']:      'CcoNH'
    }

    @staticmethod
    def getAxisDefinitionByName(name):
        for axisObj in Spectrum.ALL_AXES:
            if axisObj['label'] == name:
                return axisObj
        raise Exception("Unsupported axis definition: " + name)

    @staticmethod
    def get_cyana_type(spectrum_type):
        return Spectrum.CYANA_TYPES_DICT.get(spectrum_type)

    @staticmethod
    def getSpectrumDefinitionByName(name):
        return Spectrum.NAME_DEF_DICT.get(name, None)

    @staticmethod
    def validate_definition(definition_template):
        key = 'amino_acid_chain_shifts'
        return key in definition_template and \
               (any([i != [0] for i in definition_template[key][:-1]]) or
                min(definition_template[key][-1]) < -1 or max(definition_template[key][-1]) > 0)

    @staticmethod
    def get_peak_definition(prev_amino_acid, curr_amino_acid, definition_template):
        # later 'aa' is used as a shortcut for 'amino acid' for readability

        if Spectrum.validate_definition(definition_template):
            print("Spectrum template is probably incompatible with this peak definition generator.")

        peak_definitions = []

        shifts_prev = AminoAcid.findAminoAcidDefinition(prev_amino_acid.getAminoAcidId()) if prev_amino_acid is not None else {}
        shifts_curr = AminoAcid.findAminoAcidDefinition(curr_amino_acid.getAminoAcidId())

        main_atoms = [atom for atom in shifts_curr.keys()
                      if Atom.getAtomTypeFromLabel(atom) in definition_template['mainAtomTypes']]

        main_chain_atoms_prev = [atom for atom in shifts_prev.keys()]
        main_chain_atoms_curr = [atom for atom in shifts_curr.keys()]

        shifts = [shifts_prev, shifts_curr]
        main_chain_atoms = [main_chain_atoms_prev, main_chain_atoms_curr]

        # just as a constraint for 'for' loop not to be too long
        neighborhood = min(definition_template['neighborhood'], max(map(len, main_chain_atoms)))

        for main_atom in main_atoms:
            for main_hydrogen in shifts_curr[main_atom]:
                shifts_in_aa_chain = definition_template.get('shifts_in_amino_acid_chain', [0])
                for shift_in_aa_chain in shifts_in_aa_chain:
                    temp_shifts = shifts[shift_in_aa_chain + 1]
                    temp_main_chain_atoms = main_chain_atoms[shift_in_aa_chain + 1]
                    i = temp_main_chain_atoms.index(main_atom) if main_atom in temp_main_chain_atoms else 0
                    left, right = neighborhood, neighborhood

                    # consider HCCH-COSY when e.g. CD1 and CD2 are next to each other but both should've same neighbors
                    while (i - left >= 0 and
                           temp_main_chain_atoms[i - left][1:2] == temp_main_chain_atoms[i][1:2]) or \
                          (i - left - 1 >= 0 and
                           temp_main_chain_atoms[i - left][1:2] == temp_main_chain_atoms[i - left - 1][1:2]):
                        left += 1
                    while (i + right < len(temp_main_chain_atoms) and
                           temp_main_chain_atoms[i + right][1:2] == temp_main_chain_atoms[i][1:2]) or \
                          (i + right + 1 < len(temp_main_chain_atoms) and
                           temp_main_chain_atoms[i + right][1:2] == temp_main_chain_atoms[i + right + 1][1:2]):
                        right += 1
                    # -------------------

                    for j in range(-left, right + 1):
                        if 0 <= i + j < len(temp_main_chain_atoms):
                            neighboring_hydrogens = temp_shifts[temp_main_chain_atoms[i + j]] if Atom.getAtomTypeFromLabel(temp_main_chain_atoms[i + j]) in definition_template['hydroBondAtomTypes'] else []
                            for hydrogen in neighboring_hydrogens:
                                spectrum_axes = [axis['atom'] for axis in definition_template['axes']]
                                hydro_axes = sum([axis == 'H' for axis in spectrum_axes])
                                non_hydro_axes = definition_template['nDim'] - hydro_axes
                                axes = []
                                if non_hydro_axes >= 1:
                                    axes.append(main_atom)
                                if hydro_axes >= 1:
                                    axes.append(main_hydrogen)
                                if hydro_axes >= 2:
                                    axes.append(hydrogen)
                                if non_hydro_axes >= 2:
                                    axes.append(temp_main_chain_atoms[i + j])

                                peak_definitions.append(tuple(zip([0] * (definition_template['nDim'] - 1) + [shift_in_aa_chain], axes)))
        peak_definitions = list(map(list, list(set(peak_definitions))))
        peak_definitions.sort()

        return peak_definitions

    #
    # Constructor
    def __init__(self):
        super(Spectrum, self).__init__()
        self.axes = []
        self.order = [0, 1, 2]
        self.org_order = [0, 1, 2]
        self.mean = -1
        self.std = -1
        self.max = -1
        self.min = -1
        self.mad = -1
        self.mean_abs = -1
        self.std_abs = -1
        self.log_mean = -1
        self.max_absolute_value = -1
        self.median_abs = -1
        self.spectrum_type = None
        self.spectrum_name = None
        self.extrema = None
        self.use_true_peaks_only = False
        self.shifts = None

        self.eq_hist_data = None

    def add_config_info(self, config):
        self.spectrum_type = config.get_type()
        self.use_true_peaks_only = config.is_use_true_peaks_only()
        self.shifts = config.get_shifts()

    def update_statistics(self):

        if self.data is not None:
            abs_data = np.abs(self.data)
            self.mean = np.mean(self.data)
            self.std = np.std(self.data)
            self.max = np.max(self.data)
            self.min = np.min(self.data)
            self.mean_abs = np.mean(abs_data)
            self.std_abs = np.std(abs_data)
            self.max_absolute_value = max(abs(self.max), abs(self.min))
            self.median_abs = np.median(abs_data)
            self.median = np.median(self.data)
            self.mad = np.mean(abs_data - self.median_abs)
            self.mean_abslute_deviation = np.mean(np.abs(self.data - self.median))
            log_data = abs_data
            log_data[log_data < 1] = 1
            log_data = np.log(log_data)
            self.log_mean = np.mean(log_data)
            # self.eq_hist_data = exposure.equalize_hist(log_data)
            # if len(self.eq_hist_data.shape) == 2:
            #     self.eq_hist_data = np.asarray([self.eq_hist_data])

    #
    # Axis properties
    def get_size(self, axis_id=None):
        return self.data.shape[axis_id] if axis_id is not None else self.data.shape

    def get_num_axis(self, label):
        return sum([self.get_nucleus(x) == label for x in range(self.get_dimensionality())])

    def get_dimensionality(self):
        return len(self.axes)

    def is_2d(self):
        return self.get_nucleus(0) == Spectrum.DUMMY_AXIS['label']

    def get_num_layers(self):
        return self.data.shape[0]

    def get_axis_labels(self):
        return [self.get_nucleus(x) for x in range(3)]

    def get_axis_ppm_limits(self):
        return [axis['PPMRange'] for axis in self.axes]

    def get_nucleus(self, axis_id):
        return self.axes[axis_id]['label']

    def get_axes(self):
        return self.axes

    def get_axis_atom(self, axis_id):
        label = self.get_nucleus(axis_id)
        if label == 'NA':
            return Atom.ATOM_UNKNOWN
        elif label in ['N', 'NH', '15N']:
            return Atom.NITROGEN
        elif label in ['C', 'CH', '13C']:
            return Atom.CARBON
        elif label in ['H', 'HN', 'HC', '1H']:
            return Atom.HYDROGEN
        else:
            return Atom.ATOM_UNKNOWN

    #
    # Projections and visualizations
    def makeProjection(self, axis):
        return np.sum(self.data, axis=axis)

    def visualize_layer(self, layer_id, list_of_peaks=None, first_contour=-1):

        if first_contour == -1:
            first_contour = self.median_abs

        if self.get_dimensionality() == 2:
            if list_of_peaks is not None:
                coordinates = list_of_peaks.get_list()
                if len(coordinates.shape) < 2:
                    print("WARNING: PDF FILE NOT GENERATED " + str(coordinates.shape))
                    return

                if list_of_peaks.get_dimensionality() == 2:
                    return self.visualize_data(self.data, points_x=coordinates[:, 1], points_y=coordinates[:, 0], first_contour=first_contour)
                else:
                    return self.visualize_data(self.data, points_x=coordinates[:, 2], points_y=coordinates[:, 1], first_contour=first_contour)
            else:
                return self.visualize_data(self.data, points_x=[], points_y=[], first_contour=first_contour)

        else:
            if list_of_peaks is not None:
                coordinates = list_of_peaks.get_list_as_array()
                coordinates = coordinates[coordinates[:, 0] == layer_id, :]
                return self.visualize_data(np.squeeze(self.data[layer_id, :, :]),
                                    points_x=coordinates[:, 2], points_y=coordinates[:, 1], first_contour=first_contour)
            else:
                return self.visualize_data(np.squeeze(self.data[layer_id, :, :]), first_contour=first_contour)

    def visualize_projection(self, axis, output_path=None, list_of_peaks=None, list_of_peaks_2=None, list_of_peaks_3=None):
        data = self.makeProjection(axis)
        data /= np.max(data)

        if list_of_peaks is not None:
            coordinates = np.asarray(list_of_peaks.get_as_array()[1:], dtype=np.float64)
            print(coordinates.shape)
            coordinates = np.delete(coordinates, axis, 1)
            coordinates_2 = np.delete(np.asarray(list_of_peaks_2.get_as_array()[1:], dtype=np.float64), axis, 1) if list_of_peaks_2 is not None else []
            coordinates_3 = np.delete(np.asarray(list_of_peaks_3.get_as_array()[1:], dtype=np.float64), axis, 1) if list_of_peaks_3 is not None else []
            # if list_of_peaks_2 is None:
            #     self.visualize_data(data, points_x=coordinates[:, 1], points_y=coordinates[:, 0])
            # else:
            #     coordinates_2 = np.asarray(list_of_peaks_2.get_as_array()[1:], dtype=np.float64)
            #     coordinates_2 = np.delete(coordinates_2, axis, 1)
            self.visualize_data(data, output_path, points_x=coordinates[:, 1], points_y=coordinates[:, 0],
                                points_x_2=coordinates_2[:, 1], points_y_2=coordinates_2[:, 0],
                                points_x_3=coordinates_3[:, 1], points_y_3=coordinates_3[:, 0], number_of_contours=20)
        else:
            self.visualize_data(data, output_path=None, points_x=[], points_y=[])

    #
    # Export functions
    def save_as_ucsf(self, output_path):

        # obs = transmitter MHZ
        # sw = spectrum width

        out_axes = [a for a in self.axes if a != Spectrum.DUMMY_AXIS]
        dictionary = {'ndim': len(out_axes)}

        for axis_id in range(len(out_axes)):

            axis_def = out_axes[axis_id]
            min_ppm = axis_def['PPMRange'][0]
            max_ppm = axis_def['PPMRange'][1]
            w_freq = axis_def['frequency']
            mean = (max_ppm + min_ppm) / 2.0
            range_hz = (max_ppm - min_ppm) * w_freq

            dictionary[axis_id] = {'obs': w_freq, 'freq': True, 'label': axis_def['label'], 'complex': False,
                                   'car': mean * w_freq, 'time': False, 'size': axis_def['resolution'],
                                   'sw': range_hz, 'encoding': 'states'}

        # save spectrum
        dic = ng.sparky.create_dic(dictionary)
        ng.sparky.write(output_path, dic, np.squeeze(self.data).astype('float32'), overwrite=True)

    def add_dummy_axis(self):
        self.axes = [Spectrum.DUMMY_AXIS] + self.axes
        self.data = np.asarray([self.data])

    #
    # Other functions

    def ppm2id_scalar(self, axis_id, ppm_coord):

        axis_def = self.axes[axis_id]
        min_ppm, max_ppm = axis_def['PPMRange']
        return (max_ppm - ppm_coord) / (max_ppm - min_ppm) * axis_def['resolution']

    def ppm2id(self, ppm_coords):
        id_coords = []
        for axis_id in range(self.get_dimensionality()):
            axis_def = self.axes[axis_id]
            min_ppm, max_ppm = axis_def['PPMRange']
            if min_ppm == max_ppm:
                id_coords.append(0)
            else:
                id_coords.append((max_ppm - ppm_coords[axis_id]) / (max_ppm - min_ppm) * axis_def['resolution'])

        return id_coords

    def id2ppm(self, id_coords):
        ppm_coords = []
        for axis_id in range(self.get_dimensionality()):
            axis_def = self.axes[axis_id]
            min_ppm, max_ppm = axis_def['PPMRange']
            if min_ppm == max_ppm:
                ppm_coords.append(0)
            else:
                ppm_coords.append(max_ppm - id_coords[axis_id] * (max_ppm - min_ppm) / float(axis_def['resolution']))

        return ppm_coords

    def shift(self, ppm_coords):
        shifts = ([] if not self.is_2d() or len(ppm_coords) == 2 else [0.]) + self.shifts
        return [a - b for a, b in zip(ppm_coords, shifts)]

    def get_grid_maximum(self, id_coords):

        sign = self.data[tuple(id_coords)]
        moved = True

        while moved:

            moved = False

            for axis_id in reversed(range(self.get_dimensionality())):

                moved_dim = not (id_coords[axis_id] == 0 or id_coords[axis_id] == self.axes[axis_id]['resolution'] - 1)

                while moved_dim:
                    id_coords = tuple(id_coords)
                    points = [self.data[id_coords[:axis_id] + (id_coords[axis_id] - 1,) + id_coords[axis_id + 1:]],
                              self.data[id_coords[:axis_id] + (id_coords[axis_id],) + id_coords[axis_id + 1:]],
                              self.data[id_coords[:axis_id] + (id_coords[axis_id] + 1,) + id_coords[axis_id + 1:]]]
                    id_coords = list(id_coords)

                    if sign > 0:
                        if points[0] > points[1] and points[0] > points[2]:
                            id_coords[axis_id] -= 1
                        elif points[1] < points[2]:
                            id_coords[axis_id] += 1
                        else:
                            moved_dim = False
                    else:
                        if points[0] < points[1] and points[0] < points[2]:
                            id_coords[axis_id] -= 1
                        elif points[1] > points[2]:
                            id_coords[axis_id] += 1
                        else:
                            moved_dim = False

                    moved = moved or moved_dim

                    if id_coords[axis_id] == 0 or id_coords[axis_id] == self.axes[axis_id]['resolution'] - 1:
                        moved_dim = False

        return id_coords

    def get_centered_peak(self, id_coords):
        id_coords = tuple(id_coords)
        centered_coords = []
        for axis_id in range(self.get_dimensionality()):
            resolution = self.axes[axis_id]['resolution']
            if id_coords[axis_id] == 0 or id_coords[axis_id] == resolution - 1:
                centered_coords.append(id_coords[axis_id])
            else:
                coord = id_coords[axis_id]
                coord_minus = (coord - 1)  # % resolution
                coord_plus = (coord + 1)  # % resolution

                # def get_data(tup):
                #     Spectrum.get_mirrored_element(tup[0], tup[1], tup[2], self.data)

                x = [coord - 1, coord, coord + 1]
                # y = [get_data(id_coords[:axis_id] + (coord_minus,) + id_coords[axis_id + 1:]),
                #      get_data(id_coords[:axis_id] + (coord,) + id_coords[axis_id + 1:]),
                #      get_data(id_coords[:axis_id] + (coord_plus,) + id_coords[axis_id + 1:])]
                y = [self.data[id_coords[:axis_id] + (coord_minus,) + id_coords[axis_id + 1:]],
                     self.data[id_coords[:axis_id] + (coord,) + id_coords[axis_id + 1:]],
                     self.data[id_coords[:axis_id] + (coord_plus,) + id_coords[axis_id + 1:]]]
                coeffs = np.polyfit(x, y, deg=2)
                centered_coords.append(-coeffs[1] / (2 * coeffs[0]) if abs(coeffs[0]) > 0.0001 else coord)

        return centered_coords

    #
    # Peak filters
    def discardPeaksOnEdges(self, listOfPeaks, size=5):
        modifier = 1 if self.get_dimensionality() == 3 else 0
        listOfPeaks = listOfPeaks.get_list_excluding_range(modifier, [0, size])
        listOfPeaks = listOfPeaks.get_list_excluding_range(modifier, [self.getSize(modifier) - size, self.getSize(modifier)])
        listOfPeaks = listOfPeaks.get_list_excluding_range(modifier + 1, [0, size])
        listOfPeaks = listOfPeaks.get_list_excluding_range(modifier + 1, [self.getSize(modifier + 1) - size, self.getSize(modifier + 1)])


    @staticmethod
    def get_type_characteristics(spectrum_type):
        spectrum_def = Spectrum.NAME_DEF_DICT.get(spectrum_type)
        if spectrum_def:
            return spectrum_def['residueType'], spectrum_def['experimentClass']
        else:
            return None, None

    @staticmethod
    def parent_type_of(spectrum_type):
        types = [Spectrum.BACKBONE, Spectrum.SIDE_CHAIN, Spectrum.NOESY]
        res, exp = Spectrum.get_type_characteristics(spectrum_type)
        return types[(res == Spectrum.TYPE_SIDE_CHAIN) * (1 + (exp == Spectrum.TYPE_THROUGH_SPACE))]

    #
    # Spectrum type
    def is_noesy(self):
        _, exp = Spectrum.get_type_characteristics(self.spectrum_type)
        return exp == Spectrum.TYPE_THROUGH_SPACE

    def is_backbone(self):
        res, _ = Spectrum.get_type_characteristics(self.spectrum_type)
        return res == Spectrum.TYPE_BACKBONE

    def is_side_chain(self):
        res, exp = Spectrum.get_type_characteristics(self.spectrum_type)
        return res == Spectrum.TYPE_SIDE_CHAIN and exp != Spectrum.TYPE_THROUGH_SPACE

    def peak_list_to_grid_maximum(self, peak_list):

        result = []
        for elem in peak_list.list:
            result += [ self.get_grid_maximum([int(e) for e in elem[0:-1]]) + [elem[-1]] ]

        return PeakList(axes_names=peak_list.axes, items=result)

    def peak_list_ppm2id(self, peak_list):

        result = PeakList(peak_list.axes, items=[])

        for peak_coordinates in peak_list.get_list():
            peak_coordinates, response = peak_coordinates[:-1], peak_coordinates[-1]
            result.add_peak(self.ppm2id(peak_coordinates), response)

        return result

    def peak_list_id2ppm(self, peak_list):

        result = PeakList(peak_list.axes, items=[])

        for peak_coordinates in peak_list.get_list():
            peak_coordinates, response = peak_coordinates[:-1], peak_coordinates[-1]
            result.add_peak(self.id2ppm(peak_coordinates), response)

        return result

    def get_bundle_parallel(self, peak_list, dx, dy, dz, dxo, dyo, dzo, normalization_mode):
        num_peaks = len(peak_list)

        print('get_bundle_parallel: Processing {} peaks... '.format(num_peaks), end='')
        
        global _PARALLEL_CONTEXT
        _PARALLEL_CONTEXT = (self, dx, dy, dz, dxo, dyo, dzo, normalization_mode)
        
        start_time = timer()
        
        pool = Pool(processes=Settings.MAX_NUM_OF_CORES)
        bundles = pool.map(_parallel_work_impl, peak_list.get_list(), chunksize=max(1, num_peaks / Settings.MAX_NUM_OF_CORES))
        pool.close()
        pool.join()
        
        print('Extracted in {:.3f}s... '.format(timer() - start_time), end='')
        start_time = timer()

        bundle_data = np.array(bundles)

        labels = np.asarray([[1. - elem, elem] for elem in peak_list.get_responses()])
        print('Formatted into bundle in {:.3f}s (type: {}).'.format(timer() - start_time, bundle_data.dtype))

        peak_list_copy = PeakList.get_concatenation(peak_list)
        peak_list_copy.convert_peaks(self.id2ppm)
        peak_list_copy.reset_axes_order()

        return Bundle(data=bundle_data, labels=labels, source_labels=self.spectrum_name, peak_list=peak_list_copy, peak_list_id=peak_list)

    # @profile
    def get_extrema(self, num_peaks=-1, parallel=False):
        if self.extrema is None:
            num_layers = self.data.shape[0]

            if parallel and num_layers != 1:
                global _PARALLEL_CONTEXT
                _PARALLEL_CONTEXT = self

                pool = Pool(processes=Settings.MAX_NUM_OF_CORES)
                output_lists = pool.map(_extrema_worker, range(num_layers),
                                        chunksize=max(1, num_layers / Settings.MAX_NUM_OF_CORES))
                pool.close()
                pool.join()
            else:
                output_lists = [self.get_extrema_on_layer(x) for x in range(num_layers)]

            output_list = PeakList.get_concatenation(*output_lists)
            output_list.sort_by('Volume')
            self.extrema = output_list

        if num_peaks == -1:
            num_peaks = len(self.extrema)
        return self.extrema.get_sublist(0, num_peaks)

    # Get extrema in the image
    # @profile
    def get_extrema_on_layer(self, layer_id):

        coordinates = []

        y_max = self.data.shape[1]
        z_max = self.data.shape[2]

        for x in xrange(layer_id, layer_id+1):

            # Maxima
            tmp1 = np.logical_and(self.data[x, 2:y_max, 1:(z_max -1)] < self.data[x, 1:(y_max-1), 1:(z_max-1)],
                                  self.data[x, 0:(y_max-2), 1:(z_max -1)] < self.data[x, 1:(y_max-1), 1:(z_max-1)])
            tmp2 = np.logical_and(self.data[x, 1:(y_max-1), 2:z_max] < self.data[x, 1:(y_max-1), 1:(z_max-1)],
                                  self.data[x, 1:(y_max-1), 0: (z_max-2)] < self.data[x, 1:(y_max-1), 1:(z_max-1)])
            layer = np.logical_and(tmp1, tmp2)

            (ys,zs) = np.nonzero(layer)
            for id in xrange(0, len(ys)):
                y = ys[id]
                z = zs[id]
                if self.get_size(0) == 1 \
                    or self.get_mirrored_element(x, y+1, z+1, self.data) > self.get_mirrored_element(x+1, y+1, z+1, self.data) \
                    and self.get_mirrored_element(x, y+1, z+1, self.data) > self.get_mirrored_element(x-1, y+1, z+1, self.data):
                    coordinates.append([x, y+1, z+1, abs(self.data[x, y+1, z+1])])

            # Minima
            tmp1 = np.logical_and(self.data[x, 2:y_max, 1:(z_max -1)] > self.data[x, 1:(y_max-1), 1:(z_max-1)],
                                  self.data[x, 0:(y_max-2), 1:(z_max -1)] > self.data[x, 1:(y_max-1), 1:(z_max-1)])
            tmp2 = np.logical_and(self.data[x, 1:(y_max-1), 2:z_max] > self.data[x, 1:(y_max-1), 1:(z_max-1)],
                                  self.data[x, 1:(y_max-1), 0:(z_max-2)] > self.data[x, 1:(y_max-1), 1:(z_max-1)])
            layer = np.logical_and(tmp1, tmp2)

            (ys,zs) = np.nonzero(layer)
            for id in xrange(0, len(ys)):
                y = ys[id]
                z = zs[id]
                if self.get_size(0) == 1 \
                    or self.get_mirrored_element(x, y+1, z+1, self.data) < self.get_mirrored_element(x+1, y+1, z+1, self.data) \
                    and self.get_mirrored_element(x, y+1, z+1, self.data) < self.get_mirrored_element(x-1, y+1, z+1, self.data):
                    coordinates.append([x, y+1, z+1, abs(self.data[x, y+1, z+1])])

        return PeakList(self.get_axis_labels() + ['Volume'], items=coordinates, default_response=0.,
                        original_axes_order=self.org_order + [3])

    def get_num_extrema(self, peak_list):
        return sum([self.is_extremum(elem) for elem in peak_list])

    # Removes from PeakList negative peaks
    def filter_out_negative_peaks(self, peak_list):
        peaks = [peak for peak in peak_list if self.data[peak[0], peak[1], peak[2]] > 0]
        return PeakList(axes_names=peak_list.axes, items=peaks, original_axes_order=peak_list.org_order)

    #
    def get_volume(self, peak_cooordinates):
        suma = 0
        for x in range(int(peak_cooordinates[0])-1, int(peak_cooordinates[0])+2):
            for y in range(int(peak_cooordinates[1])-1, int(peak_cooordinates[1])+2):
                for z in range(int(peak_cooordinates[2])-1, int(peak_cooordinates[2])+2):
                    suma += abs(self.get_mirrored_element(x, y, z, self.data))
        return suma

    #
    def get_list_volume(self, peak_list):

        output = np.ones(peak_list.get_number_of_peaks())
        for i in range(0, peak_list.get_number_of_peaks()):
            output[i] = self.get_volume(peak_list.list[i][0:-1])

        return output



    # Removes from PeakList elements that do not lie in extrema
    def filter_out_nonextrema(self, peak_list):
        peaks = [peak for peak in peak_list if self.is_extremum(peak)]
        return PeakList(axes_names=peak_list.axes, items=peaks, original_axes_order=peak_list.org_order)

    # Removes from PeakList elements that lie on the diagonal
    def filter_out_diagonal_peaks(self, peak_list):
        axis_ratio = self.get_size(1) / self.get_size(2)
        peaks = [peak for peak in peak_list if abs(peak[1] - axis_ratio * peak[2]) > 2]
        return PeakList(axes_names=peak_list.axes, items=peaks, original_axes_order=peak_list.org_order)

    # Removes from PeakList elements that lie in the water region
    def filter_out_peaks_in_water_range(self, peak_list, horizontal=True, vertical=True):
        peaks = peak_list.list
        for axis_id in range(self.get_dimensionality()):

            range_h = self.get_water_range_id(axis_id)

            if not range_h == [] and (horizontal and axis_id == 1 or vertical and axis_id == 2):
                peaks = [peak for peak in peaks if not range_h[0] <= peak[axis_id] <= range_h[1]]

        return PeakList(axes_names=peak_list.axes, items=peaks, original_axes_order=peak_list.org_order)

    def pull_peaks_to_limit(self, peak_list):
        for axis_id, limit in enumerate(self.get_size()):
            for peak in peak_list:
                peak[axis_id] = min(peak[axis_id], limit - 1)

    def get_water_range_id(self, axis_id):

        if self.get_axis_atom(axis_id) == Atom.HYDROGEN:
            return [self.ppm2id_scalar(axis_id, 5.8), self.ppm2id_scalar(axis_id, 3.8)]
        else:
            return []

    def count_positive_peaks(self, peak_list):
        return sum( 1 for elem in peak_list.list if self.data[elem[0], elem[1], elem[2]] > 0)

    def is_extremum(self, coords):
        # TODO:
        # use ">=" instead of ">"?

        elem = self.get_mirrored_element(coords[0], coords[1], coords[2], self.data)
        isMaximum = \
            (self.get_size(0) == 1 or elem > self.get_mirrored_element(coords[0] + 1, coords[1], coords[2], self.data) and elem > self.get_mirrored_element(coords[0] - 1, coords[1], coords[2], self.data)) \
            and elem > self.get_mirrored_element(coords[0], coords[1] + 1, coords[2], self.data) \
            and elem > self.get_mirrored_element(coords[0], coords[1] - 1, coords[2], self.data) \
            and elem > self.get_mirrored_element(coords[0], coords[1], coords[2] + 1, self.data) \
            and elem > self.get_mirrored_element(coords[0], coords[1], coords[2] - 1, self.data)
        isMinimum = \
            (self.get_size(0) == 1 or elem < self.get_mirrored_element(coords[0] + 1, coords[1], coords[2], self.data) and elem < self.get_mirrored_element(coords[0] - 1, coords[1], coords[2], self.data)) \
            and elem < self.get_mirrored_element(coords[0], coords[1] + 1, coords[2], self.data) \
            and elem < self.get_mirrored_element(coords[0], coords[1] - 1, coords[2], self.data) \
            and elem < self.get_mirrored_element(coords[0], coords[1], coords[2] + 1, self.data) \
            and elem < self.get_mirrored_element(coords[0], coords[1], coords[2] - 1, self.data)
        return isMinimum or isMaximum

    #
    # Allows to access single element of an array. It index is outside array range, a mirroring will be used
    # to provide final value.
    @staticmethod
    def get_mirrored_element(x, y, z, data):

        if x < 0:
            x = min(data.shape[0]-1, abs(x))
        elif x >= data.shape[0]:
            x = max(0, 2 * (data.shape[0] - 1) - x)
        
        if y < 0:
            y = min(data.shape[1]-1, abs(y))
        elif y >= data.shape[1]:
            y = max(0, 2 * (data.shape[1] - 1) - y)
            
        if z < 0:
            z = min(data.shape[2]-1, abs(z))
        elif z >= data.shape[2]:
            z = max(0, 2 * (data.shape[2] - 1) - z)

        return data[x, y, z]

    def rescale_spectrum(self, newSize):

        currSize = np.asarray(self.data.shape, dtype=np.float32)
        self.data = ndimage.zoom( self.data, np.asarray(newSize, dtype=np.float32) / currSize, order=1 )
        for elem in range(0, len(newSize)):
            self.axes[elem]['resolution'] = self.data.shape[elem]

    #
    #
    def get_cube(self, x, y, z, dx, dy, dz, data=None):
        if data is None:
            data = self.data

        shape = data.shape
        if dx <= x < shape[0] - dx and dy <= y < shape[1] - dy and dz <= z < shape[2] - dz:
            x, y, z = int(x), int(y), int(z)
            return np.array(data[x - dx:x + dx + 1, y - dy:y + dy + 1, z - dz:z + dz + 1])

        result = np.zeros([2 * dx + 1, 2 * dy + 1, 2 * dz + 1], dtype=np.float32)

        for xs in xrange(0, result.shape[0]):
            x_tmp = x + xs - dx
            for ys in xrange(0, result.shape[1]):
                y_tmp = y + ys - dy
                for zs in xrange(0, result.shape[2]):
                    result[xs, ys, zs] = Spectrum.get_mirrored_element(x_tmp, y_tmp, z + zs - dz, data)

        return result

