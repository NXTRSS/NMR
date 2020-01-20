import csv
import random
import warnings
from itertools import compress
import scipy.io as sio
from nmr_commons.sequence.Atom import Atom
import numpy as np
import re
from copy import deepcopy


# CLASS FUNCTIONS AND PROPERTIES:
#   - Stores 2D or 3D list of peaks
#   - Basic list modifications and filtering
#   - IO operations
#   - Data analysis
# Last line-by-line check: 20-10-2016
# Tests: nmr_commons/tests/PeakListTest.py


class PeakList:
    DUMMY_AXIS_NAME = 'NA'

    ####
    # Constructor
    def __init__(self, axes_names, items=None, default_response=1.0, original_axes_order=None):
        self.list = deepcopy(items) if items else []
        self.axes = list(axes_names)
        self.org_order = list(range(self.get_dimensionality()) if original_axes_order is None else original_axes_order)

        for peak in self.list:
            if len(peak) == self.get_dimensionality():
                peak.append(default_response)

    ####
    # Peak addition, removal, modification
    def add_peak(self, peak, response=1.0):
        self.list.append([float(i) for i in peak] + [response])

    def add_peaks(self, *peak_lists):
        for peak_list in peak_lists:
            if self.axes == peak_list.axes or not all(self.axes):
                self.list.extend(deepcopy(peak_list.get_list()))
            else:
                warnings.warn('Peak axes mismatch when calling add_peaks. Some peaks are not added to the list.')

    def convert_peaks(self, coord_conversion_function):
        self.list = [coord_conversion_function(peak[:-1]) + [peak[-1]] for peak in self.list]

    def permute(self, permutation):
        self.list = np.asarray(self.list)[permutation].tolist()

    def remove_duplicates(self):
        self.list = sorted(map(list, set(map(tuple, self.list))))

    def set_axes(self, axes):
        self.axes = deepcopy(axes)

    def set_responses(self, responses):
        if isinstance(responses, (int, float, np.int, np.float)):
            responses = [responses] * len(self)
        if len(responses) != len(self):
            warnings.warn('Dimension mismatch in set_responses.')
        else:
            for peak, response in zip(self.list, responses):
                peak[-1] = response

    def shift_to_ppm_limits(self, limits, epsilon=0):
        if len(limits) != self.get_dimensionality():
            warnings.warn('Dimension mismatch in shift_to_ppm_limits.')
        else:
            peaks_changed = 0
            for peak in self.list:
                peak_cp = list(peak)
                for i in range(self.get_dimensionality()):
                    min_ppm, max_ppm = limits[i]
                    while not min_ppm - epsilon <= peak[i] <= max_ppm + epsilon:
                        peak[i] += (max_ppm - min_ppm) * np.sign(max_ppm - peak[i])
                peaks_changed += peak != peak_cp
            return peaks_changed

    def unfold(self, area, limits, shift_axis, multiplicity=1, eps=1e-4):
        if len(limits) != self.get_dimensionality():
            warnings.warn('Dimension mismatch in shift_beyond_ppm_limits.')
        else:
            peaks_changed = 0
            min_ppm, max_ppm = limits[shift_axis]
            for peak in self.list:
                if all(area[axis][0] - eps <= peak[axis] <= area[axis][1] + eps for axis in range(self.get_dimensionality())):
                    peak[shift_axis] += multiplicity * (max_ppm - min_ppm)
                    peaks_changed += 1
            return peaks_changed

    def unfold_based_on_ideal(self, limits, ideal_list, unfold_axis):
        ideal_list = ideal_list.get_list_in_range_for_unfolding(limits, unfold_axis)
        areas = ideal_list.get_unfold_areas(limits, unfold_axis)
        peaks_changed = 0
        for area in areas:
            peaks_changed += self.unfold(area['limits'], limits, unfold_axis,
                                         area['shift_multiplicity'])
        return peaks_changed

    def shift_by(self, offset):
        offset += [0]  # no offset for response, not using 'append' as it would modify method argument
        self.list = [[a + b for a, b in zip(peak, offset)] for peak in self.list]

    def shift_by_for_unfolded(self, offset, limits, unfold_axis, eps=0.001):
        min_ppm, max_ppm = limits[unfold_axis]
        spectrum_width = max_ppm - min_ppm
        for index, peak in enumerate(self.list):
            shifted = False
            for i in range(len(offset)):
                if min_ppm + i*spectrum_width - eps < peak[unfold_axis] < min_ppm + (i+1)*spectrum_width + eps:
                    if not shifted:
                        temp_offset = offset[i] + [0]
                        self.list[index] = [a + b for a, b in zip(peak, temp_offset)]
                        shifted = True

    def shuffle_peaks(self):
        random.shuffle(self.list)

    def sort_by(self, column, ascending=False):
        if isinstance(column, str):
            column = (self.axes + ['response']).index(column)
        self.list.sort(key=lambda peak: (1 if ascending else -1) * peak[column])

    ####
    # Column manipulation
    def add_column(self, value, header):
        if not isinstance(value, list):
            value = [value] * len(self)
        if len(value) == len(self):
            self.axes.append(header)
            for peak, new_value in zip(self.list, value):
                peak.append(new_value)

    def adjust_axes_order(self):
        weighted_axes = sorted([(-Atom.get_atom_weight(axis), i, axis) for i, axis in
                                zip(range(self.get_dimensionality()), self.axes)])
        _, axis_order, ordered_axes = zip(*weighted_axes)

        for i in range(self.get_dimensionality()):
            self.org_order[axis_order[i]] = i

        self.reorder_axes(axis_order)

    def insert_dummy_axis(self):
        if self.axes[0] != PeakList.DUMMY_AXIS_NAME:
            self.org_order.insert(0, -1)
            self.org_order = [order + 1 for order in self.org_order]
            self.axes.insert(0, PeakList.DUMMY_AXIS_NAME)
            for peak in self.list:
                peak.insert(0, 0.)

    def remove_column(self, axis):
        if isinstance(axis, str):
            axis = self.axes.index(axis)
        if axis < 0:
            axis += self.get_dimensionality()
        if axis >= self.get_dimensionality():
            warnings.warn('Axis index out of range. Column not removed.')
        for peak in self.list:
            del peak[axis]
        del self.axes[axis]

    def remove_dummy_axis(self):
        if self.axes[0] == PeakList.DUMMY_AXIS_NAME:
            self.remove_column(0)
            self.org_order = self.org_order[1:]

    def reorder_axes(self, new_order):
        if all([type(axis) == str for axis in new_order]):
            warnings.warn('Reordering using string labels is unstable when labels are non-unique.')
            if sorted(new_order) != sorted(self.axes):
                atoms_axes = [Atom.get_universal_atom_label(x) for x in self.axes]
                if sorted(new_order) != sorted(atoms_axes):
                    return
                for i in range(self.get_dimensionality()):
                    self.swap_axes(i, atoms_axes.index(new_order[i]))
                    atoms_axes = [Atom.get_universal_atom_label(x) for x in self.axes]
            else:
                for i in range(self.get_dimensionality()):
                    self.swap_axes(i, self.axes.index(new_order[i]))
        else:
            if len(new_order) == self.get_dimensionality() \
                    and min(new_order) == 0 and max(new_order) == self.get_dimensionality() - 1:
                self.axes = [self.axes[i] for i in new_order]
                self.list = [[peak[i] for i in new_order] + [peak[-1]] for peak in self.list]

    def reset_axes_order(self):
        self.reorder_axes(self.org_order)
        self.org_order = list(range(self.get_dimensionality()))

    def swap_axes(self, axis1, axis2):
        for peak in self.list:
            peak[axis1], peak[axis2] = peak[axis2], peak[axis1]
        self.axes[axis1], self.axes[axis2] = self.axes[axis2], self.axes[axis1]

    ####
    # Getters for modified lists
    def get_binarized_list(self, threshold=0):
        items = [peak[:-1] + [1 if peak[-1] >= threshold else 0] for peak in self.list]
        return PeakList(axes_names=self.axes, items=items, original_axes_order=self.org_order)

    @staticmethod
    def get_concatenation(peak_list, *peak_lists):
        p = PeakList(axes_names=peak_list.axes, original_axes_order=peak_list.org_order)
        p.add_peaks(peak_list, *peak_lists)
        return p

    def get_filtered_list(self, threshold, rule='above'):
        return PeakList(axes_names=self.axes,
                        items=[elem for elem in self.list if ( elem[-1] >= threshold and rule == 'above' or elem[-1] < threshold and rule == 'below')],
                        original_axes_order=self.org_order)

    def get_list_excluding_peaks(self, peak_list):
        responses = set(self.get_responses())
        if len(responses) == 1:
            result_list = self._perform_fast_excluding(peak_list, list(responses)[0])
        else:
            result_list = [peak for peak in self.list if not peak_list.is_in_list(peak)]
        return PeakList(axes_names=self.axes, items=result_list, original_axes_order=self.org_order)

    def _perform_fast_excluding(self, peak_list, response):
        peaks = {tuple(x[:-1]) for x in self}
        to_exclude = {tuple(x[:-1]) for x in peak_list}
        return [list(peak) + [response] for peak in peaks.difference(to_exclude)]

    def get_list_excluding_range(self, axis_id, excluded_range):
        modifier = 1 if self.is_dummy_axis_present() else 0
        items = [peak for peak in self.list if not excluded_range[0] <= peak[axis_id + modifier] <= excluded_range[1]]
        return PeakList(axes_names=self.axes, items=items, original_axes_order=self.org_order)

    def get_lists_grouped_by(self, axis_id):
        pass

    def get_sublist(self, first, last=None):
        if isinstance(first, (list, tuple, np.ndarray)) and last is None:
            if all([isinstance(x, bool) or x.dtype == np.bool for x in first]):
                items = list(compress(self.list, first))
            elif all([isinstance(x, int) or x.dtype == np.int for x in first]):
                items = [self.list[ind] for ind in first]
            else:
                raise Exception('Invalid argument "first" in get_sublist. '
                                'The list should contain only bool or only int values.')
        else:
            items = self.list[first:last]
        return PeakList(axes_names=self.axes, items=items, original_axes_order=self.org_order)

    def get_list_in_range_for_unfolding(self, limits, unfold_axis, eps=0.001):
        ideal_list = self
        for axis in range(len(self.axes)):
            if axis != unfold_axis:
                min_ppm, max_ppm = limits[axis]
                ideal_list = ideal_list.get_list_excluding_range(axis, [min_ppm-100, min_ppm+eps])
                ideal_list = ideal_list.get_list_excluding_range(axis, [max_ppm-eps, max_ppm+100])
        return ideal_list

    ####
    # Getters
    def get_column(self, axis):
        if isinstance(axis, str):
            axis = self.axes.index(axis)
        return [peak[axis] for peak in self.list]

    def get_list(self):
        return self.list

    def get_as_array(self):
        return np.concatenate(([self.axes + ['response']], self.list))

    def get_dimensionality(self):
        return len(self.axes)

    def get_number_of_peaks(self):
        warnings.warn('peak_list.get_number_of_peaks() is deprecated. Use len(peak_list) instead.')
        return len(self.list)

    def __len__(self):
        return len(self.list)

    def __iter__(self):
        return iter(self.list)

    def __getitem__(self, item):
        return self.list[item]

    def get_responses(self):
        return self.get_column(-1)

    def get_unfold_areas(self, limits, shift_axis):
        max_coord = np.max(self.get_column(shift_axis))
        min_ppm, max_ppm = limits[shift_axis]
        spectrum_width = max_ppm - min_ppm
        areas = []
        shift_multiplicity = 0
        lower_bound = max_ppm
        while lower_bound < max_coord:
            shift_multiplicity += 1
            upper_bound = lower_bound + spectrum_width
            in_range = [peak for peak in self.list if lower_bound < peak[shift_axis] < upper_bound]
            area_limits = [[np.min([peak[axis] for peak in in_range]), np.max([peak[axis] for peak in in_range])] for
                           axis in range(self.get_dimensionality())]
            area_limits[shift_axis][0] -= (shift_multiplicity * spectrum_width)
            area_limits[shift_axis][1] -= (shift_multiplicity * spectrum_width)
            area = {'limits': area_limits, 'shift_multiplicity': shift_multiplicity}
            areas += [area]
            lower_bound = upper_bound
        return areas

    def get_unfold_axis(self, limits):
        unfold_axis = -1
        max_peaks_beyond_limits = -1
        for axis in range(len(self.axes)):
            min_ppm, max_ppm = limits[axis]
            peaks_beyond_limits = sum(1 for peak in self.list if peak[axis] < min_ppm or peak[axis] > max_ppm)
            if peaks_beyond_limits > max_peaks_beyond_limits:
                unfold_axis = axis
                max_peaks_beyond_limits = peaks_beyond_limits
        return unfold_axis

    def round_coordinates(self, precision=0):
        self.convert_peaks(lambda peak: [round(coord, precision) for coord in peak])
        if precision == 0:
            self.convert_peaks(lambda peak: [int(round(coord)) for coord in peak])  # stupid python2 needs it

    ####
    # IO operations
    def save_as_csv(self, path, response=True):
        with open(path, 'w') as f:
            writer = csv.writer(f, delimiter=',')
            writer.writerow(self.axes + (['response'] if response else []))
            writer.writerows(self.list if response else [peak[:-1] for peak in self.list])

    def save_as_mat(self, path):
        sio.savemat(path, {'list': self.list, 'axes': self.axes})

    def save_as_sparky(self, path, order=(0, 1, 2)):
        dim = sum([axis.startswith('w') for axis in self.axes])
        if dim == 0:
            dim = self.get_dimensionality()
        with open(path, 'w') as f:
            f.write('{0: >17}'.format('Assignment '))
            for i in range(dim):
                f.write('{0: >11}'.format('w' + str(i + 1)))
            for i in range(dim, self.get_dimensionality()):
                f.write('{0: >13}'.format(self.axes[i]))
            f.write('\n')

            for peak in self.list:
                f.write('{0: >17}'.format('?-' * (dim - 1) + '?'))
                for i in range(dim):
                    f.write('{:11.3f}'.format(peak[order[i]]))
                for i in range(dim, self.get_dimensionality()):
                    f.write('{:13}'.format(peak[i]))
                f.write('\n')

    def save_as_xeasy(self, path, header):
        with open(path, 'w') as f:
            f.write(header)

            for i, peak in enumerate(self.list):
                f.write('{:>7} '.format(i))
                for j in range(self.get_dimensionality()):
                    f.write('{:7.3f} '.format(peak[j]))
                f.write('{} {}   {:.3e}  {:.3e} {} {} {:>5} {:>5}'.format(1, 'U', peak[-1], 0, 'e', 0, 0, 0))
                if self.get_dimensionality() == 3:
                    f.write(' {:>5}'.format(0))
                f.write('\n')

    @staticmethod
    def load_from_xeasy(path):
        axes_names = []
        with open(path, 'r') as f:
            next_line = f.readline()
            while next_line.startswith('#'):
                if next_line.startswith('#SPECTRUM'):
                    axes_names = next_line.split()[2:]
                next_line = f.readline()
            peak_list = PeakList(axes_names=axes_names,
                                 items=[[float(x) for x in line.split()[1:4]] for line in [next_line] + f.readlines() if
                                        line != '\n'])
            return peak_list

    @staticmethod
    def load_from_mat(path):
        mat_struct = sio.loadmat(path)
        result = PeakList([str(a).strip() for a in mat_struct['axes']])
        result.list = mat_struct['list'].tolist()
        return result

    @staticmethod
    def load_from_csv(path, includes_response='auto', default_response=1.0):
        with open(path, 'r') as f:
            reader = csv.reader(f, delimiter=',')
            header = next(reader)
            if includes_response == 'auto':
                includes_response = header[-1] == 'response'

            peak_list = PeakList(axes_names=header[:-1 if includes_response else len(header)])
            peak_list.list = [list(map(float, row)) + ([] if includes_response else [default_response])
                              for row in reader]

            # handle old NMR-Datatool format
            if peak_list.axes[0] == '2D':
                peak_list.axes[0] = PeakList.DUMMY_AXIS_NAME
            return peak_list

    @staticmethod
    def load_from_dumpling_csv(path):
        with open(path, 'r') as f:
            [f.readline() for _ in range(2)]  # two lines of header
            axis_names = f.readline()
            peak_list = PeakList(axes_names=re.findall("\d+\D", axis_names))
            f.readline()  # header

            while True:
                line = f.readline().split(';')[0].split(',')
                if line != ['']:
                    peak_list.add_peak([float(x) for x in line])
                else:
                    break

            return peak_list

    @staticmethod
    def load_from_sparky(path):
        with open(path, 'r') as f:
            header = f.readline()
            peak_list = PeakList(axes_names=header.split()[1:],
                                 items=[[float(x) for x in line.split()[1:]] for line in f.readlines() if line != '\n'])
            return peak_list

    @staticmethod
    def load_from_xpk(path):
        with open(path, 'r') as f:
            [f.readline() for _ in range(1)]
            axis_name =[axis[0:-1] for axis in f.readline().split()]
            peak_list = PeakList(axes_names=axis_name)
            [f.readline() for _ in range(4)]
            if len(axis_name) == 3:
                while True:
                    line = f.readline().split()
                    if line != []:
                        peaks = [line[2], line[9], line[16]]
                        peak_list.add_peak([float(x) for x in peaks])
                    else:
                        break
            if len(axis_name) == 2:
                while True:
                    line = f.readline().split()
                    if line != []:
                        peaks = [line[2], line[9]]
                        peak_list.add_peak([float(x) for x in peaks])
                    else:
                        break

        return peak_list

    ###
    # Data analysis


    ####
    # Others
    @staticmethod
    def from_array(array):
        return PeakList(axes_names=array[0, :-1].tolist(), items=np.asarray(array[1:], dtype=np.float64).tolist())

    @staticmethod
    def get_shift_from_file(shift_path):
        with open(shift_path, 'r') as f:
            shifts = [line.split()[1:] for line in f if line.strip()]
            shifts = [[float(shift) for shift in axis_shifts] for axis_shifts in shifts]
        shifts = map(list, zip(*shifts))
        return shifts

    def is_dummy_axis_present(self):
        return self.axes[0] == PeakList.DUMMY_AXIS_NAME

    def is_in_list(self, other_peak, only_coords=False):
        return any([peak[:-1] == other_peak[:len(other_peak) if only_coords else -1] for peak in self.list])

    def to_string(self):
        return str(self.axes) + '\n' + '\n'.join(map(str, self.list))

    def exclude_peaks_outside_spectrum_range(self, spectrum):
        r1 = spectrum.axes[0]['resolution']
        r2 = spectrum.axes[1]['resolution']
        r3 = spectrum.axes[2]['resolution']
        output = []

        for elem in self.list:

            if all( [ a[0] < a[1] for a in zip(elem[0:2], [r1, r2, r3]) ] ) \
                    and all( [ a >= 0 for a in elem[0:2] ] ):
                output.append(elem)
            else:
                print("Peak excluded: " + str(elem))

        return PeakList(self.axes, items=output)






