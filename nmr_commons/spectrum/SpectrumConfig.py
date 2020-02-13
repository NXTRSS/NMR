import ConfigParser


class SpectrumConfig:

    def __init__(self):
        self.parser = ConfigParser.ConfigParser()
        self.parser.add_section('shifts')
        self.parser.add_section('info')
        self.parser.add_section('unfold_shifts')

    def set_shifts(self, w1, w2, w3):
        self.parser.set('shifts', 'w1', w1)
        self.parser.set('shifts', 'w2', w2)
        self.parser.set('shifts', 'w3', w3)

    def get_shifts(self):
        shifts = [float(self.parser.get('shifts', 'w1')), float(self.parser.get('shifts', 'w2'))]
        if self.parser.has_option('shifts', 'w3'):
            shifts.append(float(self.parser.get('shifts', 'w3')))
        return shifts

    def get_type(self):
        return self.parser.get('info', 'type')

    def set_type(self, type):
        return self.parser.set('info', 'type', type)

    def is_use_true_peaks_only(self):
        return self.parser.has_option('info', 'true_peaks_only') and self.parser.get('info', 'true_peaks_only')

    def save(self, file_path):
        with open(file_path, 'w') as f:
            self.parser.write(f)

    def get_unfold_shifts(self, axis_name):
        return [float(x) for x in self.parser.get('unfold_shifts', axis_name).split(',')]

    def get_unfold_axis(self):
        if self.parser.has_option('unfold_axis', 'axis'):
            return self.parser.get('unfold_axis', 'axis')
        return None

    @staticmethod
    def load(file_path):
        output = SpectrumConfig()
        output.parser.read(file_path)
        return output
