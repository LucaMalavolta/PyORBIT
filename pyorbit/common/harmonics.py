from pyorbit.common.abstract_common import *

class CommonHarmonics(AbstractCommon):

    model_class = 'harmonics'

    ''' all the possible parameters that can be assigned to the harmonics are listed here'''

    parameters_dictionary = {
        'P': # Main period of the harmonic series
            {
                'bounds': [0.5, 1e03],
                'priors': ['Uniform', []],
                'spaces': 'Log_Base2',
                'fixed' : None,
                'unit': 'days',
            },
        'T0': # reference time of the harmic series
            {
                'bounds': [-1e02, 1e02],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.00,
                'unit': None,
            },
        'phase': # Phase offset, shared among all the harmonics
            {
                'bounds': [0.0, 360.],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.00,
                'unit': 'degree',
            },
    }
    recenter_pams = {'phase'}

    for i_order in range(1,6):
        # Amplitude of the sine component of the i_order harmonic
        parameters_dictionary['amp_S'+repr(i_order)] = {
                'bounds': [1e-06, 1e06],
                'priors': ['Uniform', []],
                'spaces': 'Log_Base2',
                'fixed' : 0.00,
                'unit': 'as input',
        }
        # Phase of the sine component of the i_har harmonic
        parameters_dictionary['pha_S'+repr(i_order)] = {
                'bounds': [0.0, 360. / i_order],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.00,
                'unit': 'degree',
        }
        recenter_pams.add('pha_S'+repr(i_order))

        # Amplitude of the cosine component of the i_har harmonic
        parameters_dictionary['amp_C'+repr(i_order)] = {
                'bounds': [1e-06, 1e06],
                'priors': ['Uniform', []],
                'spaces': 'Log_Base2',
                'fixed' : 0.00,
                'unit': 'as input',
        }
        # Phase of the cosine component of the i_har harmonic
        parameters_dictionary['pha_C'+repr(i_order)] = {
                'bounds': [0.0, 360. / i_order],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.00,
                'unit': 'degree',
        }
        recenter_pams.add('pha_C'+repr(i_order))

    default_fixed = {}

