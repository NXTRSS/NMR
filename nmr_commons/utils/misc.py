import os


def flatten(source_list):
    return [j for i in source_list for j in i]


def get_spectra_names(main_dir):
    extension = '.ucsf'
    return [(main_dir, elem[:-len(extension)]) for elem in os.listdir(main_dir) if elem.endswith(extension)]
