from pyorbit.common.abstract_common import *

class CommonHarmonics(AbstractCommon):

    model_class = 'harmonics'

    ''' all the possible parameters that can be assigned to the harmonics are listed here'''
    list_pams = {
        'P',  # Periodicity of the first harmonics
        'T0', # Referecne time for the sinusoid
        'phase',  # Phase offset, shared amond all the harmonics
        'amp_S1',  # Amplitude of the sine component of the first harmonic
        'amp_S2',  # Amplitude of the sine component of the second harmonic
        'amp_S3',  # Amplitude of the sine component of the third harmonic
        'amp_S4',  # Amplitude of the sine component of the fourth harmonic
        'amp_S5',  # Amplitude of the sine component of the fifth harmonic
        'amp_C1',  # Amplitude of the cosine component of the first harmonic
        'amp_C2',  # Amplitude of the cosine component of the second harmonic
        'amp_C3',  # Amplitude of the cosine component of the third harmonic
        'amp_C4',  # Amplitude of the cosine component of the fourth harmonic
        'amp_c5',  # Amplitude of the cosine component of the fifth harmonic
        'pha_S1',  # Phase of the sine component of the first harmonic
        'pha_S2',  # Phase of the sine component of the second harmonic
        'pha_S3',  # Phase of the sine component of the third harmonic
        'pha_S4',  # Phase of the sine component of the fourth harmonic
        'pha_S5',  # Phase of the sine component of the fifth harmonic
        'pha_C1',  # Phase of the cosine component of the first harmonic
        'pha_C2',  # Phase of the cosine component of the second harmonic
        'pha_C3',  # Phase of the cosine component of the third harmonic
        'pha_C4',  # Phase of the cosine component of the fourth harmonic
        'pha_c5',  # Phase of the cosine component of the fifth harmonic
    }

    """These default boundaries are used when the user does not define them in the yaml file"""
    default_bounds = {
        'P': [0.5, 1000.0],
        'T0': [-10 ** 2, 10 ** 2],  # Amplitude of the sine component of the first harmonic
        'phase': [0.0, 2 * np.pi],  # Phase offset, shared amond all the harmonics
        'amp_S1': [10 ** -6, 10 ** 6],  # Amplitude of the sine component of the first harmonic
        'amp_S2': [10 ** -6, 10 ** 6],  # Amplitude of the sine component of the second harmonic
        'amp_S3': [10 ** -6, 10 ** 6],  # Amplitude of the sine component of the third harmonic
        'amp_S4': [10 ** -6, 10 ** 6],  # Amplitude of the sine component of the fourth harmonic
        'amp_S5': [10 ** -6, 10 ** 6],  # Amplitude of the sine component of the fifth harmonic
        'amp_C1': [10 ** -6, 10 ** 6],  # Amplitude of the cosine component of the first harmonic
        'amp_C2': [10 ** -6, 10 ** 6],  # Amplitude of the cosine component of the second harmonic
        'amp_C3': [10 ** -6, 10 ** 6],  # Amplitude of the cosine component of the third harmonic
        'amp_C4': [10 ** -6, 10 ** 6],  # Amplitude of the cosine component of the fourth harmonic
        'amp_c5': [10 ** -6, 10 ** 6],  # Amplitude of the cosine component of the fifth harmonic
        'pha_S1': [0.0, 2 * np.pi],
        'pha_S2': [0.0, 2 * np.pi / 2.],
        'pha_S3': [0.0, 2 * np.pi / 3.],
        'pha_S4': [0.0, 2 * np.pi / 4.],
        'pha_S5': [0.0, 2 * np.pi / 5.],
        'pha_C1': [0.0, 2 * np.pi],
        'pha_C2': [0.0, 2 * np.pi / 2.],
        'pha_C3': [0.0, 2 * np.pi / 3.],
        'pha_C4': [0.0, 2 * np.pi / 4.],
        'pha_C5': [0.0, 2 * np.pi / 5.],
    }

    """ These default priors are used when the user does not define them in the yaml file
        The boundaries are not included in this dictionary, because it is likely that the user will specify his
        preferred boundaries without changing the priors
    """
    default_priors = {
        'P': ['Uniform', []],
        'T0': ['Uniform', []],
        'phase': ['Uniform', []],
        'amp_S1': ['Uniform', []],
        'amp_S2': ['Uniform', []],
        'amp_S3': ['Uniform', []],
        'amp_S4': ['Uniform', []],
        'amp_S5': ['Uniform', []],
        'amp_C1': ['Uniform', []],
        'amp_C2': ['Uniform', []],
        'amp_C3': ['Uniform', []],
        'amp_C4': ['Uniform', []],
        'amp_C5': ['Uniform', []],
        'pha_S1': ['Uniform', []],
        'pha_S2': ['Uniform', []],
        'pha_S3': ['Uniform', []],
        'pha_S4': ['Uniform', []],
        'pha_S5': ['Uniform', []],
        'pha_C1': ['Uniform', []],
        'pha_C2': ['Uniform', []],
        'pha_C3': ['Uniform', []],
        'pha_C4': ['Uniform', []],
        'pha_C5': ['Uniform', []],
    }

    default_spaces = {
        'P': 'Log_Base2',
        'T0': 'Linear',
        'phase': 'Linear',
        'amp_S1': 'Log_Base2',
        'amp_S2': 'Log_Base2',
        'amp_S3': 'Log_Base2',
        'amp_S4': 'Log_Base2',
        'amp_S5': 'Log_Base2',
        'amp_C1': 'Log_Base2',
        'amp_C2': 'Log_Base2',
        'amp_C3': 'Log_Base2',
        'amp_C4': 'Log_Base2',
        'amp_C5': 'Log_Base2',
        'pha_S1': 'Linear',
        'pha_S2': 'Linear',
        'pha_S3': 'Linear',
        'pha_S4': 'Linear',
        'pha_S5': 'Linear',
        'pha_C1': 'Linear',
        'pha_C2': 'Linear',
        'pha_C3': 'Linear',
        'pha_C4': 'Linear',
        'pha_C5': 'Linear',
    }

    default_fixed = {}

    recenter_pams = {'phase'}
