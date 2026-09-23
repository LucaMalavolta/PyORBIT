from pyorbit.common.abstract_common import *
from pyorbit.keywords_definitions import *

class CommonSpectrograph(AbstractCommon):
    model_class = 'spectrograph'
    ''' all the possible parameters that can be assigned to the spectrograph are listed here'''

    parameters_dictionary = {
    'line_broadening':
        {
            'bounds': [0.0, 5.0],
            'priors': ['Uniform', []],
            'spaces': 'Linear',
            'fixed' : 1.5,
            'unit': 'km/s',
        },
    'line_contrast':
        {
            'bounds': [0.0, 1.0],
            'priors': ['Uniform', []],
            'spaces': 'Linear',
            'fixed' : 0.5,
            'unit': 'relative depth',
        },
    'instrumental_broadening': # FWHM of the spectra line
        {
            'bounds': [0.0, 10.0],
            'priors': ['Uniform', []],
            'spaces': 'Linear',
            'fixed' : 1.108,
            'unit': 'km/s',
        },
    'ccf_broadening':
        {
            'bounds': [0.0, 10.0],
            'priors': ['Uniform', []],
            'spaces': 'Linear',
            'fixed' : 1.5,
            'unit': 'km/s',
        },
    'measured_ccf_width': # FWHM
        {
            'bounds': [0.0, 200.0],
            'priors': ['Uniform', []],
            'spaces': 'Linear',
            'fixed' : 1.6,
            'unit': 'km/s',
        },
    }

    recenter_pams = {}

    def __init__(self, *args, **kwargs):
        super(CommonSpectrograph, self).__init__(*args, **kwargs)

        self.rv_min = -20.00 # km/s
        self.rv_max = 20.00 # km/s
        self.rv_step = 0.5 # km/s
        self.use_stellar_lines = True

    def initialize_model(self, mc, **kwargs):

        for keyword in keywords_rv_min:
            self.rv_min = kwargs.get(keyword, self.rv_min)
        for keyword in keywords_rv_max:
            self.rv_max = kwargs.get(keyword, self.rv_max)
        for keyword in keywords_rv_step:
            self.rv_step = kwargs.get(keyword, self.rv_step)

        for keyword in keywords_use_stellar_lines:
            self.use_stellar_lines = kwargs.get(keyword, self.use_stellar_lines)

    def print_info(self):
        print("*** spectrograph {0:s} global parameters:".format(self.common_ref))

        print("    Use stellar line parameters (in common for all instruments): ",self.use_stellar_lines)
        print("    CCF RV starting value: {0:f.3}  km/s".format(self.rv_min))
        print("    CCF RV end value:      {0:f.3}  km/s".format(self.rv_max))
        print("    CCF RV step:      {0:f.3}  km/s".format(self.step))
        print("    Remember to put a prior or fix the instrumental_broadening")