from abstract_common import *


class CommonActivity(AbstractCommon):

    model_class = 'activity'

    ''' all the possible parameters that can be assigned to a planet are listed here'''
    list_pams = {
        'Prot': 'U',  # Rotational period of the star
        'Pdec': 'U',  # Decay timescale of activity
        'Oamp': 'LU',  # Granulation of activity
        'Hamp': 'U',  # Amplitude of the signal in the covariance matrix
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
        'cel_a': [0.00000001, 1000000.0],
        'cel_b': [0.00000001, 1000000.0],
        'cel_c': [0.00000001, 1000000.0]
    }

    recenter_pams = {}
