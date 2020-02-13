import matplotlib.pyplot as plt
import numpy as np

from nmr_commons.peak_picking.PeakList import PeakList
from nmr_commons.spectrum.Data import Data
import matplotlib
import _pickle as pickle
import scipy.signal as signal
import scipy.ndimage.filters as filters


class Bundle(Data):

    # Constructor
    def __init__(self, data=None, labels=None, source_labels=None, peak_list=None, peak_list_id=None):
        super(Bundle, self).__init__()
        self.data = data
        self.labels = None
        self.source_labels = None
        self.spectrum_axes_dict = {}
        self.peak_list = peak_list
        self.peak_list_id = peak_list_id

        if labels is not None:
            self.set_labels(labels)
        if source_labels is not None:
            self.set_source_labels(source_labels)
        if self.peak_list is None:
            self.peak_list = PeakList(axes_names=[None, None, None])
        elif all(self.peak_list.axes):
            if self.source_labels is not None:
                for spectrum in set(self.source_labels):
                    if spectrum not in self.spectrum_axes_dict:
                        self.spectrum_axes_dict[spectrum] = list(self.peak_list.axes)

    def __len__(self):
        return len(self.data) if self.data is not None else None

    def set_labels(self, value):
        if isinstance(value, np.ndarray):
            self.labels = value
        elif isinstance(value, list):
            self.labels = np.asarray(value)
        else:
            self.labels = np.ones(self.data.shape[0]) * value

    def set_source_labels(self, value):
        if isinstance(value, (np.ndarray, list)):
            self.source_labels = value
        else:
            self.source_labels = np.array(self.data.shape[0] * [value])

    def visualize_bundle(self, subplot_rows=3, subplot_cols=3, save_path=None, additional_labels=None, title_fig='', limit=-1):

        if limit == -1:
            limit = self.data.shape[0]
        else:
            limit = min(limit, self.data.shape[0])

        number_of_subplots = subplot_rows * subplot_cols

        matplotlib.rc('xtick', labelsize=6)
        matplotlib.rc('ytick', labelsize=6)

        for elem in xrange(0, limit, number_of_subplots):

            f = plt.figure(facecolor='white', figsize=(8.27, 11.69))

            for i in xrange(0, number_of_subplots):
                global_index = elem + i
                a = f.add_subplot(subplot_rows, subplot_cols, i+1)
                if self.labels.shape[0] == self.data.shape[0]:
                    if global_index < limit:
                        self.visualize_data(
                            np.squeeze(self.data[global_index, :, :, 0]),
                            title=("Example {} : {} {:.2f}".format(
                                global_index,
                                self.labels[global_index],
                                additional_labels[global_index] if additional_labels is not None else '_')
                            ),
                            axes_obj=a)
                else:
                    self.visualize_data(np.squeeze(self.data[(global_index):, :, 0]), axes_obj=a)

            plt.suptitle(title_fig)

            if save_path is not None:
                f.savefig(save_path + str(elem+i) + ".png", bbox_inches='tight')

        return f

    def add_bundle(self, bundle):
        if self.data is None:
            self.data = np.array(bundle.data)
            self.labels = np.array(bundle.labels)
            self.source_labels = np.array(bundle.source_labels)
        else:
            self.data = np.concatenate((self.data, bundle.data), axis=0)
            self.labels = np.concatenate((self.labels, bundle.labels), axis=0)
            self.source_labels = np.concatenate((self.source_labels, bundle.source_labels), axis=0)
            if self.data.shape[0] != len(self.labels):
                raise Exception("Number of labels in the bundle does not match number of observations")

        self.peak_list.add_peaks(bundle.peak_list)
        self.spectrum_axes_dict.update(bundle.spectrum_axes_dict)

    def save(self, path):
        bundle_dict = [[key] + values for key, values in self.spectrum_axes_dict.iteritems()]
        np.savez_compressed(path, data=self.data, labels=self.labels, source_labels=self.source_labels,
                            peak_list=self.peak_list.get_as_array(), bundle_dict=bundle_dict)

    @staticmethod
    def load(path):
        arrays = np.load(path)
        bundle = Bundle(data=arrays['data'], labels=arrays['labels'], source_labels=arrays['source_labels'],
                        peak_list=PeakList.from_array(arrays['peak_list']))
        bundle.spectrum_axes_dict = dict([(row[0], row[1:].tolist()) for row in arrays['bundle_dict']])
        return bundle

    def permute(self):
        permutation = np.random.permutation(len(self.labels))
        data, labels, source_labels = self.data[permutation], self.labels[permutation], self.source_labels[permutation]
        del [self.data, self.labels, self.source_labels]
        self.data, self.labels, self.source_labels = data, labels, source_labels
        self.peak_list.permute(permutation)

    def get_subbundle_indices(self, idx):
        bundle = Bundle(data=self.data[idx,:,:,:], peak_list=self.peak_list.get_sublist(idx))
        bundle.set_labels([self.labels[id] for id in idx])
        bundle.set_source_labels([self.source_labels[id] for id in idx])
        return bundle



    def get_subbundle(self, first, last):
        bundle = Bundle(data=self.data[first:last], peak_list=self.peak_list.get_sublist(first, last))
        bundle.set_labels(self.labels[first:last])
        bundle.set_source_labels(self.source_labels[first:last])
        spectra = set(bundle.source_labels)
        for key, value in self.spectrum_axes_dict.iteritems():
            if key in spectra:
                bundle.spectrum_axes_dict[key] = list(value)
        return bundle

    def get_subbundle_ignore(self, first, minibatch_size, ignored_spectra):
        bundle = Bundle(data=np.zeros((minibatch_size, self.data.shape[1], self.data.shape[2], self.data.shape[3])),
                        labels=np.zeros((minibatch_size, 2)), source_labels=['']*minibatch_size)
        entry = 0
        for elem_id in xrange(first, len(self)):

            if self.source_labels[elem_id] not in ignored_spectra:

                bundle.data[entry] = self.data[elem_id, :, :, :]
                bundle.labels[entry] = self.labels[elem_id]
                bundle.source_labels[entry] = self.source_labels[elem_id]
                peak = self.peak_list[elem_id]
                bundle.peak_list.add_peak(peak[:-1], peak[-1])

                curr_spectrum = bundle.source_labels[entry]
                if curr_spectrum not in bundle.spectrum_axes_dict:
                    bundle.spectrum_axes_dict[curr_spectrum] = list(self.spectrum_axes_dict[curr_spectrum])

                entry += 1

            if entry == minibatch_size:
                return bundle, elem_id

        return bundle, len(self)


    def convolve(self, gauss_size):
        for i in range(0, self.data.shape[0]):
            self.data[i, :, :, 0] = filters.gaussian_filter(np.squeeze(self.data[i, :, :, :]), sigma=gauss_size)

    def get_minibatches(self, batch_size, ignored=None):

        if ignored is None:
            for batch_index in range(0, len(self), batch_size):
                yield self.get_subbundle(batch_index, batch_index + batch_size)
        else:
            first = 0
            while first < len(self):
                bundle, index = self.get_subbundle_ignore(first, batch_size, ignored)
                first = index
                yield bundle
