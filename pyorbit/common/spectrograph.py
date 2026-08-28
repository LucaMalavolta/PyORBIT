from pyorbit.common.abstract_common import *

class CommonSpectrograph(AbstractCommon):
    model_class = 'spectrograph'
    ''' all the possible parameters that can be assigned to the spectrograph are listed here'''

    parameters_dictionary = {
    'natural_broadening':
        {
            'bounds': [0.0, 5.0],
            'priors': ['Uniform', []],
            'spaces': 'Linear',
            'fixed' : 1.5,
            'unit': 'km/s',
        },
    'natural_contrast':
        {
            'bounds': [0.0, 1.0],
            'priors': ['Uniform', []],
            'spaces': 'Linear',
            'fixed' : 0.5,
            'unit': 'relative depth',
        },
    'instrumental_broadening':
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
