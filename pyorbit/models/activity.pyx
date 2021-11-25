from pyorbit.models.abstract_common import *


class CommonActivity(AbstractCommon):

    model_class = 'activity'

    ''' all the possible parameters that can be assigned to a planet are listed here'''
    list_pams = {
        # quasi-periodic 
        'Prot',  # Rotational period of the star
        'Pdec',  # Decay timescale of activity
        'Oamp',  # Granulation of activity
        'Hamp',  # Amplitude of the signal in the covariance matrix
        'Camp',  # Amplitude of the derivative of the kernel
        #celerite 2 parameters
        'grn_period', # undamped period of the granulation
        'grn_sigma', # the standard deviation of the process
        'rot_sigma',
        'rot_fmix',
        'rot_log10_Q0',
        'rot_log10_deltaQ',
        # GP framework parameter
        'Vc',  # GP framework parameter
        'Vr',  # GP framework parameter
        'Lc',  # GP framework parameter
        'Bc',  # GP framework parameter
        'Br',  # GP framework parameter
         # Matern kernel
        'matern32_sigma', # Matern kernel,  sigma parameter
        'matern32_rho', # Matern kernel,  rho parameter
        'matern32_log10_sigma', # Matern kernel,  sigma parameter
        'matern32_log10_rho', # Matern kernel,  rho parameter
        # sinusoid-like activity
        'P',  # Period
        'K',  # Sinusoid semi-amplitude
        'f',  # Sinusoid curve phase
        # legacy parameters for the second implementation of celerite
        'Hamp_factor',
        'Q0',  # celerite rotation term parameter
        'deltaQ', # celerite rotation term parameter
        'ln_Q0',  # celerite rotation term parameter
        'ln_deltaQ', # celerite rotation term parameter
        'mix', # celerite rotation term parameter
        'amp', # celerite rotation term parameter
        # legacy parameters for the first implementation of celerite
        'cel_factor',
        'cel_a',  # celerite term A
        'cel_b',  # celerite term B
        'cel_c',  # celerite term C
    }

    """These default boundaries are used when the user does not define them in the yaml file"""
    default_bounds = {
        'Prot': [1.0, 1000.0],
        'Pdec': [1.0, 10000.0],
        'Oamp': [0.0001, 2.0],
        'Hamp': [0.00000001, 1000000.0],
        'Camp': [0.00000001, 1000000.0],
        #
        'grn_period': [1.0, 1000.0],
        'grn_sigma': [0.00000001, 1000000.0],
        'rot_sigma': [0.00000001, 1000000.0],
        'rot_fmix': [0.001, 1.000],
        'rot_Q0': [0.001, 1000.000],
        'rot_deltaQ': [0.001, 1000.000],
        'rot_log10_Q0':[-10., 10.0],
        'rot_log10_deltaQ': [-10., 10.0],
        #
        'Vc': [-500.0, 500.0],
        'Vr': [-500.0, 500.0],
        'Lc': [-500.0, 500.0],
        'Bc': [-500.0, 500.0],
        'Br': [-500.0, 500.0],
        #
        'P': [0.4, 100000.0],
        'K': [0.5, 2000.0],
        'f': [0.0, 2 * np.pi],
        #
        'matern32_sigma': [0.00000001, 1000000.0],
        'matern32_rho': [0.001, 1000.00],
        #
        'matern32_log10_sigma': [-6.0, 6.0],
        'matern32_log10_rho': [-3.0, 3.0],
        #
        'Hamp_factor': [0.01, 10.0],
        'Q0': [0.001, 1000.000],
        'deltaQ': [0.001, 1000.000],
        'ln_Q0': [-10., 10.0],
        'ln_deltaQ': [-10., 10.0],
        'mix': [0.001, 1000.000],
        'amp': [0.0001, 1000.0],
        #
        'cel_factor': [0.00000001, 1000000.0],
        'cel_B': [0.00000001, 1000000.0],
        'cel_C': [0.00000001, 1000000.0],
    }

    """ These default priors are used when the user does not define them in the yaml file
        The boundaries are not included in this dictionary, because it is likely that the user will specify his
        preferred boundaries without changing the priors
    """
    default_priors = {
        'Prot': ['Uniform', []],
        'Pdec': ['Uniform', []],
        'Oamp': ['Uniform', []],
        'Hamp': ['Uniform', []],
        'Camp': ['Uniform', []],
        #
        'grn_period': ['Uniform', []], # undamped period of the granulation
        'grn_sigma': ['Uniform', []], # the standard deviation of the process
        'rot_sigma': ['Uniform', []],
        'rot_fmix': ['Uniform', []],
        'rot_Q0': ['Uniform', []],
        'rot_deltaQ': ['Uniform', []],
        'rot_log10_Q0': ['Uniform', []],
        'rot_log10_deltaQ': ['Uniform', []],
        #
        'Vc': ['Uniform', []],
        'Vr': ['Uniform', []],
        'Lc': ['Uniform', []],
        'Bc': ['Uniform', []],
        'Br': ['Uniform', []],
        #
        'matern32_sigma': ['Uniform', []],
        'matern32_rho': ['Uniform', []],
        'matern32_log10_sigma': ['Uniform', []],
        'matern32_log10_rho': ['Uniform', []],
        #
        'P': ['Uniform', []],
        'K': ['Uniform', []],
        'f': ['Uniform', []],
        #
        'Hamp_factor': ['Uniform', []],
        'Q0': ['Uniform', []],
        'deltaQ': ['Uniform', []],
        'ln_Q0': ['Uniform', []],
        'ln_deltaQ': ['Uniform', []],
        'mix': ['Uniform', []],
        'amp': ['Uniform', []],
        #
        'cel_factor': ['Uniform', []],
        'cel_B': ['Uniform', []],
        'cel_C': ['Uniform', []],
    }

    default_spaces = {
        'Prot': 'Linear',
        'Pdec': 'Linear',
        'Oamp': 'Logarithmic',
        'Hamp': 'Linear',
        'Camp': 'Linear',
        #
        'grn_period': 'Linear',
        'grn_sigma': 'Logarithmic',
        'rot_sigma': 'Logarithmic',
        'rot_fmix': 'Linear',
        'rot_Q0': 'Logarithmic',
        'rot_deltaQ': 'Logarithmic',
        'rot_log10_Q0': 'Linear',
        'rot_log10_deltaQ': 'Linear',
        #
        'Vc': 'Linear',
        'Vr': 'Linear',
        'Lc': 'Linear',
        'Bc': 'Linear',
        'Br': 'Linear',
        #
        'matern32_sigma':'Linear',
        'matern32_rho': 'Linear',
        #
        'matern32_log10_sigma':'Linear',
        'matern32_log10_rho': 'Linear',
        #
        'P': 'Logarithmic',
        'K': 'Logarithmic',
        'f': 'Linear',
        #
        'Hamp_factor': 'Linear',
        'Q0': 'Logarithmic',
        'deltaQ': 'Logarithmic',
        'ln_Q0': 'Linear',
        'ln_deltaQ': 'Linear',
        'mix': 'Logarithmic',
        'amp': 'Logarithmic',
        #
        'cel_factor': 'Linear',
        'cel_B': 'Logarithmic',
        'cel_C': 'Logarithmic',
    }

    default_fixed = {}

    recenter_pams = {'f'}
