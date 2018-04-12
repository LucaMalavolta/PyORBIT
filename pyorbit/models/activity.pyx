from abstract_common import *


class CommonActivity(AbstractCommon):

    model_class = 'activity'

    ''' all the possible parameters that can be assigned to a planet are listed here'''
    list_pams = {
        'Prot': 'U',  # Rotational period of the star
        'Pdec': 'U',  # Decay timescale of activity
        'Oamp': 'LU',  # Granulation of activity
        'Hamp': 'U',  # Amplitude of the signal in the covariance matrix
        'P': 'LU',  # Period, log-uniform prior
        'K': 'LU',  # RV semi-amplitude, log-uniform prior
        'f': 'U',  # RV curve phase, log-uniform
        'cel_a': 'LU',  # celerite term A
        'cel_b': 'LU',  # celerite term B
        'cel_c': 'LU',  # celerite term C
    }

    """These default boundaries are used when the user does not define them in the yaml file"""
    default_bounds = {
        'Prot': [1.0, 1000.0],
        'Pdec': [1.0, 10000.0],
        'Oamp': [0.0001, 2.0],
        'Hamp': [0.00000001, 1000000.0],
        'P': [0.4, 100000.0],
        'K': [0.5, 2000.0],
        'f': [0.0, 2 * np.pi],
        'cel_a': [0.00000001, 1000000.0],
        'cel_b': [0.00000001, 1000000.0],
        'cel_c': [0.00000001, 1000000.0]
    }

    """ These default priors are used when the user does not define them in the yaml file
        The boundaries are not included in this dictionary, because it is likely that the user will specify his
        preferred boundaries without changing the priors
    """
    default_priors = {
        'Prot': ['Uniform', []],
        'Pdec': ['Uniform', []],
        'Oamp': ['Jeffreys', []],
        'Hamp': ['Uniform', []],
        'P': ['Jeffreys', []],
        'K': ['ModifiedJeffreys', [1.0]],
        'f': ['Uniform', []],
        'cel_a': ['Jeffreys', []],
        'cel_b': ['Jeffreys', []],
        'cel_c': ['Jeffreys', []]
    }

    recenter_pams = {'f'}
