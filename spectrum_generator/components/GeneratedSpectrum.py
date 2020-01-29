import scipy.ndimage as nd
from timeit import default_timer as timer
from scipy.ndimage.filters import gaussian_filter
from nmr_commons.spectrum.Spectrum import Spectrum
from nmr_commons.peak_picking.PeakList import PeakList
import matplotlib.pyplot as plt
from spectrum_generator.components.rbm.lasagne_rbm import processImage
from spectrum_generator.components.rbm.lasagne_rbm import processEntireImage
import spectrum_generator.components.rbm
import random
from nmr_environment.Settings import Settings
import numpy as np
import copy


class GeneratedSpectrum(Spectrum):

    SIGMA = 1
    SIZE = 15
    CENTER = int(SIZE/2.0 + 0.5)
    GAUSS_3D = np.zeros([SIZE, SIZE, SIZE])
    GAUSS_3D[CENTER, CENTER, CENTER] = 1
    GAUSS_3D = gaussian_filter(GAUSS_3D, SIGMA)
    GAUSS_3D = GAUSS_3D/GAUSS_3D[CENTER, CENTER, CENTER]
    GAUSS_2D = np.squeeze(GAUSS_3D[CENTER, :, :])

    def __init__(self, spectrum_def):
        super(GeneratedSpectrum, self).__init__()

        if isinstance(spectrum_def, (tuple, list)):
            self.axes = None
            self.spectrum_def = None
            self.resolution = list(spectrum_def)
        else:
            self.axes = copy.deepcopy(spectrum_def['axes'])
            self.spectrum_def = copy.deepcopy(spectrum_def)
            self.resolution = [ax['resolution'] for ax in self.axes]

        self.data = np.zeros(self.resolution)

    @staticmethod
    def from_sparky_spectrum(spectrum):
        """
        Creates GeneratedSpectrum object using config provided in real SparkySpectrum object.
        """
        gen_spectrum = GeneratedSpectrum(spectrum.get_size())
        gen_spectrum.axes = copy.deepcopy(spectrum.axes)
        return gen_spectrum

    @staticmethod
    def get_default_generator_parameters():

        # Generator parameters
        generator_parameters = {}
        generator_parameters['mean_intensity'] = 2577909.03812
        generator_parameters['std_intensity'] = 400000
        # generator_parameters['peaks_seed_type'] = 'synthetic'
        # generator_parameters['peaks_seed'] = 'gaussian'
        # generator_parameters['artifacts_seed_type'] = 'synthetic'
        # generator_parameters['artifacts_seed'] = 'spike_and_slab'
        generator_parameters['peaks_seed_type'] = 'file'
        generator_parameters['peaks_seed'] = 'seed_2'
        generator_parameters['artifacts_seed_type'] = 'file'
        generator_parameters['artifacts_seed'] = 'seed_noise'
        generator_parameters['mean_intensity_mul'] = 0.7
        generator_parameters['std_intensity_mul'] = 0.001

        return generator_parameters

    #
    #
    def insert_image(self, center_position, i):

        image_size = np.asarray(i.shape)

        image_left_down_corner = np.ceil(np.floor(center_position) - image_size / 2.0)
        image_right_upper_corner = np.ceil(np.floor(center_position) + image_size / 2.0)
        image_left_down_corner_truncated = [int(x) for x in np.maximum([0] * self.get_dimensionality(), image_left_down_corner)]
        image_right_upper_corner_truncated = [int(x) for x in np.minimum(self.data.shape, image_right_upper_corner)]

        local_left_down_corner = [int(x) for x in np.abs(image_left_down_corner_truncated - image_left_down_corner)]
        local_right_upper_corner = [int(x) for x in image_size - np.abs(image_right_upper_corner_truncated - image_right_upper_corner)]

        if any(np.asarray(image_right_upper_corner_truncated) <= np.asarray(image_left_down_corner_truncated)):
            print("Warning: image can not be places at position " + str(center_position))
            return


        image_data = self.data[image_left_down_corner_truncated[0]:image_right_upper_corner_truncated[0],
                    image_left_down_corner_truncated[1]:image_right_upper_corner_truncated[1],
                    image_left_down_corner_truncated[2]:image_right_upper_corner_truncated[2]]
        new_data = i[local_left_down_corner[0]:local_right_upper_corner[0],
                    local_left_down_corner[1]:local_right_upper_corner[1],
                    local_left_down_corner[2]:local_right_upper_corner[2]]
        image_larger = np.abs(image_data) > np.abs(new_data)

        self.data[image_left_down_corner_truncated[0]:image_right_upper_corner_truncated[0],
            image_left_down_corner_truncated[1]:image_right_upper_corner_truncated[1],
            image_left_down_corner_truncated[2]:image_right_upper_corner_truncated[2]] = image_data * image_larger + new_data * (1-image_larger)




    #
    @staticmethod
    def get_noisy_gaussians():

        X_temp = np.zeros((1,1,32,32))
        m_dist = 4
        sigma_1 = 1
        sigma_2 = 1
        for i1 in range(0,32):
            for i2 in range(0,32):

                for g1 in range(4, 32, m_dist):
                    for g2 in range(4, 32, m_dist):
                        gaussValue = 1 / ( 2 * np.pi * sigma_1 * sigma_2 ) * np.exp(-1/2 * ((i1-g1)**2/(sigma_1**2)+(i2-g2)**2/(sigma_2**2)))
                        X_temp[:,:,i1,i2] = np.maximum(X_temp[:,:,i1,i2], gaussValue)

        #X_temp = X_temp + (np.random.randn(X_temp.shape[0], X_temp.shape[1], X_temp.shape[2], X_temp.shape[3])/noise if noise > 0 else 0)

        return X_temp

    # Noisy image
    @staticmethod
    def get_noisy_image(noise=500.0):
        X_temp = np.zeros((1,1,61,61))
        X_temp = X_temp + (np.random.randn(X_temp.shape[0], X_temp.shape[1], X_temp.shape[2], X_temp.shape[3])/noise if noise > 0 else 0)
        return X_temp

    # Generate realistic map
    def generate_realistic_map(self, peak_list, gener_params, generate_water=False, generate_noise=False):

        t = timer()

        # Build data structure for response map
        spectrum_resolution = self.get_size()
        print("GENERATING PEAK IMAGES ...")
        I = processEntireImage(K=len(peak_list.get_list()), gener_params=gener_params, model_mode=1,
                               model_name=Settings.RBM_PATH)
        # Generate peak images
        i = 0
        for peak_coordinates in peak_list.get_list():
            peak_coordinates, response = peak_coordinates[:-1], peak_coordinates[-1]
            self.insert_image(peak_coordinates, I[i, :, :])
            i += 1

        # Generate water images
        if generate_water:
            for peakId in range(0,60):
                peak_coordinates_X = self.data.shape[1] - np.abs(np.random.normal(0, 1))
                peak_coordinates_Y = np.random.randint(0, self.data.shape[0])

                i = processImage(GeneratedSpectrum.get_noisy_gaussian(noise=400.0, sigma_1=0.5, sigma_2=1.0, m=15.5))
                self.insert_image([peak_coordinates_Y, peak_coordinates_X], i)

                print(peak_coordinates_Y)
                print(peak_coordinates_X)
                print( self.data.shape[1] - peak_coordinates_X)

        # Generate noise
        if generate_noise:
            d_y = I.shape[2] / 2.0
            d_z = I.shape[3] / 2.0
            y = spectrum_resolution[1]
            z = spectrum_resolution[2]
            n_y = int(np.ceil(y/d_y))+1
            n_z = int(np.ceil(z/d_z))+1
            print(int(n_y*n_z))
            noise = n_y*n_z
            print("GENERATING NOISE IMAGES ...")
            I_noise = processEntireImage(K=int(noise), gener_params=gener_params, model_mode=2, model_name=Settings.RBM_PATH_N)

            idx = 0
            i = 0
            for i_y in range(0, n_y):
                for i_z in range(0, n_z):
                    idy = int((d_y-1)/2 + i_y*d_y) + np.random.uniform( - d_y/3.0, d_y/3.0)
                    idz = int((d_z-1)/2 + i_z*d_z) + np.random.uniform( - d_z/3.0, d_z/3.0)
                    self.insert_image( [idx, idy, idz], I_noise[i,:,:] )
                    i += 1

        self.update_statistics()

    #
    # Response map generation
    def generate_response_map(self, peak_list, convert_from_ppm=False, verbose=False, verbose_warn_only=False):
        t = timer()

        invisible_peaks = 0

        # Convert peak coordinates
        if convert_from_ppm:
            peak_list.convert_peaks(self.ppm2id)
        peak_list.round_coordinates()
        for peak in peak_list:
            peak, response = peak[:-1], peak[-1]
            if all(0 <= elem < max_ind for elem, max_ind in zip(peak, self.resolution)):
                self.data[tuple(peak)] = response
            else:
                invisible_peaks += 1

        if verbose:
            print('Generating response map... ({} peak(s) will not be visible)'.format(invisible_peaks))
        if verbose_warn_only and invisible_peaks > 0:
            print('WARNING: {} peak(s) will not be visible on generated response map.'.format(invisible_peaks))

        # Convolve initial response map with gaussian
        weights_scheme = GeneratedSpectrum.GAUSS_3D if self.get_dimensionality() == 3 else GeneratedSpectrum.GAUSS_2D
        data = nd.convolve(self.data, weights_scheme, origin=1.0)

        if verbose:
            print('Response map generated in {:.2f}s.'.format(timer() - t))

        self.data = data  # * (1.0/np.max(data[:]))
        self.update_statistics()
