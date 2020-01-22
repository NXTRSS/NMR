import os

class FileUtils:

    @staticmethod
    def get_trajectory_names(path, extension):
        return [(path, elem[:-len(extension)]) for elem in os.listdir(path) if elem.endswith(extension)]

    @staticmethod
    def create_directory(path):
        if not os.path.exists(path):
            os.makedirs(path)

    @staticmethod
    def get_spectra_names(main_dir):
        extension = '.ucsf'
        return [(main_dir, elem[:-len(extension)]) for elem in os.listdir(main_dir) if elem.endswith(extension)]