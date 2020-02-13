from nmr_commons.spectrum.SpectrumReader import SpectrumReader
from nmr_commons.spectrum.Spectrum import Spectrum
from nmr_commons.peak_picking.PeakList import PeakList
from nmr_commons.sequence.Atom import Atom
import nmrglue as ng
import numpy as np
import os
from nmr_commons.spectrum.SpectrumConfig import SpectrumConfig


class SparkySpectrum(Spectrum):

    ####
    # Constructor
    def __init__(self, file_path, config_path=None, low_mem=False, extrema_path=None, adjust_axis_order=False):
        super(SparkySpectrum, self).__init__()

        self.dic, self.data = ng.sparky.read(file_path) #SpectrumReader.read(file_path, low_mem)
        self.file_path = file_path
        self.update_statistics()

        # spectrum name
        head, tail = os.path.split(file_path)
        self.spectrum_name = tail[0:-5]

        # extrema path
        if extrema_path is not None:
            self.extrema = PeakList.load_from_mat(extrema_path)

        # spectrum config info
        if config_path is not None:
            config = SpectrumConfig.load(config_path)
            self.add_config_info(config)

        # Build axis
        num_dim = self.dic['naxis']

        for w in range(num_dim):
            uc = ng.sparky.make_uc(self.dic, self.data, dim=w)
            ppm_range = uc.ppm_limits()
            new_axis = {'label': self.dic['w{}'.format(w + 1)]['nucleus'],
                        'PPMRange': [ppm_range[1], ppm_range[0]],
                        'resolution': self.data.shape[w],
                        'frequency': self.dic['w{}'.format(w + 1)]['spectrometer_freq']}
            self.axes.append(new_axis)

        if num_dim == 2:
            self.add_dummy_axis()

        if adjust_axis_order:
            self.adjust_axes_order()

    def add_dummy_axis(self):
        super(SparkySpectrum, self).add_dummy_axis()
        self.dic['w3'] = self.dic['w2']
        self.dic['w2'] = self.dic['w1']
        self.dic['w1'] = None

    def to_string(self):
        return "Spectrum " + self.spectrum_name + " of type " + self.spectrum_type

    def adjust_axes_order(self):
        if self.extrema is not None:
            self.extrema.adjust_axes_order()

        weighted_axes = sorted([(-Atom.get_atom_weight(axis), i, axis) for i, axis in
                                zip(range(self.get_dimensionality()), self.get_axis_labels())])
        _, axis_order, ordered_axes = zip(*weighted_axes)

        self.reorder_axes(axis_order)

    def reorder_axes(self, new_order):
        if len(new_order) == self.get_dimensionality() \
                and min(new_order) == 0 and max(new_order) == self.get_dimensionality() - 1:
            self.axes = [self.axes[i] for i in new_order]
            self.data = self.data.transpose(new_order)
            dic_copy = {'w1': self.dic['w1'], 'w2': self.dic['w2'], 'w3': self.dic['w3']}
            for i in range(self.get_dimensionality()):
                self.dic['w{}'.format(i + 1)] = dic_copy['w{}'.format(new_order[i] + 1)]

            self.order = tuple([self.order[i] for i in new_order])
            org_order_copy = list(self.org_order)
            for i in range(self.get_dimensionality()):
                self.org_order[new_order[i]] = org_order_copy[i]

    def reset_axes_order(self):
        self.reorder_axes(self.org_order)
        self.org_order = list(range(self.get_dimensionality()))
        self.order = tuple(self.org_order)

    def get_extrema(self, num_peaks=-1, parallel=False):
        old_order = self.order
        self.adjust_axes_order()
        extrema = super(SparkySpectrum, self).get_extrema(num_peaks, parallel)
        self.reset_axes_order()
        self.reorder_axes(old_order)
        return extrema
