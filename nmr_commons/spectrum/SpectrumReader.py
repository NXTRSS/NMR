import os
from warnings import warn
import nmrglue as ng
import numpy as np


class SpectrumReader:

    def __init__(self, file_path):
        self.file_path = file_path
        self.dic = None
        self.data = None

    @staticmethod
    def read(file_path, low_mem):
        if low_mem:
            return ng.sparky.read_lowmem(file_path)
        else:
            reader = SpectrumReader(file_path)
            reader._read()
            return reader.dic, reader.data
        
    def _read(self):
        seek_pos = os.stat(self.file_path).st_size
        with open(self.file_path, 'rb') as f:
            self.dic = ng.sparky.fileheader2dic(ng.sparky.get_fileheader(f))

            size_mismatch = seek_pos != self.dic['seek_pos'] and self.dic['seek_pos'] != 0
            if size_mismatch:
                warn('Bad file size in header: {} (file) vs {} (header)'.format(seek_pos, self.dic['seek_pos']))

            for i in range(self.dic['naxis']):
                self.dic['w' + str(i + 1)] = ng.sparky.axisheader2dic(ng.sparky.get_axisheader(f))

            bytes_excluded = seek_pos - self.dic['seek_pos'] if size_mismatch else None
            self.data, excluded = self._get_data(f, bytes_excluded)

            if self.dic['naxis'] == 2:
                self._transform_data_2d()
            elif self.dic['naxis'] == 3:
                self._transform_data_3d()
            else:
                raise ValueError('Unknown dimensionality: {}'.format(self.dic['naxis']))

            if excluded is not None:
                warn('Excluded part:\n' + excluded)

    def _transform_data_2d(self):
        len_y = self.dic['w1']['npoints']
        len_x = self.dic['w2']['npoints']
        tile_y = self.dic['w1']['bsize']
        tile_x = self.dic['w2']['bsize']
        self.data = ng.sparky.untile_data2D(self.data, (tile_y, tile_x), (len_y, len_x))

    def _transform_data_3d(self):
        len_z = self.dic['w1']['npoints']
        len_y = self.dic['w2']['npoints']
        len_x = self.dic['w3']['npoints']
        tile_z = self.dic['w1']['bsize']
        tile_y = self.dic['w2']['bsize']
        tile_x = self.dic['w3']['bsize']
        self.data = ng.sparky.untile_data3D(self.data, (tile_z, tile_y, tile_x), (len_z, len_y, len_x))

    @staticmethod
    def _get_data(f, bytes_excluded):
        data = f.read()
        if bytes_excluded:
            include, exclude = data[:-bytes_excluded], data[-bytes_excluded:]
        else:
            include, exclude = data, None
        return np.frombuffer(include, dtype='>f4'), exclude
