from pyorbit.common.abstract_common import *
from pyorbit.keywords_definitions import *

class CommonActivity(AbstractCommon):

    model_class = 'activity'

    ''' all the possible parameters that can be assigned to activity models are listed here'''
    parameters_dictionary = {
        # quasi-periodic
        'Prot':  # Rotational period of the star
            {
                'bounds': [1.0, 1000.0],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : None,
                'unit': 'days',
            },
        'Pdec':  # Decay timescale of activity
            {
                'bounds': [1.0, 1000.0],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : None,
                'unit': 'days',
            },
        'Pcyc':  # Decay timescale of activity
            {
                'bounds': [50.0, 10000.0],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : None,
                'unit': 'days',
            },
        'Oamp':  # Granulation of activity
            {
                'bounds': [0.0001, 2.0],
                'priors': ['Uniform', []],
                'spaces': 'Log_Base2',
                'fixed' : None,
                'unit': 'as input',
            },
        'Hamp':  # Amplitude of the signal in the covariance matrix
            {
                'bounds': [0.00000001, 1000000.0],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : None,
                'unit': 'as input',
            },
        'Camp':  # Amplitude of the derivative of the kernel
            {
                'bounds': [0.00000001, 1000000.0],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : None,
                'unit': 'as input',
            },
        #celerite 2 parameters
        'sho_period': # undamped period of the granulation
            {
                'bounds': [1.0, 1000.0],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : None,
                'unit': 'days',
            },
        'sho_tau': # the standard deviation of the process
            {
                'bounds': [0.00000001, 1000000.0],
                'priors': ['Uniform', []],
                'spaces': 'Log_Base2',
                'fixed' : None,
                'unit': 'as input',
            },
        'sho_sigma': # the standard deviation of the process
            {
                'bounds': [0.00000001, 1000000.0],
                'priors': ['Uniform', []],
                'spaces': 'Log_Base2',
                'fixed' : None,
                'unit': 'as input',
            },
        'grn_period': # undamped period of the granulation
            {
                'bounds': [1.0, 1000.0],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : None,
                'unit': 'days',
            },
        'grn_sigma': # the standard deviation of the process
            {
                'bounds': [0.00000001, 1000000.0],
                'priors': ['Uniform', []],
                'spaces': 'Log_Base2',
                'fixed' : None,
                'unit': 'as input',
            },
        'rot_sigma':
            {
                'bounds': [0.00000001, 1000000.0],
                'priors': ['Uniform', []],
                'spaces': 'Log_Base2',
                'fixed' : None,
                'unit': 'as input',
            },
        'rot_fmix':
            {
                'bounds': [0.001, 1.0],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : None,
                'unit': 'as input',
            },
        'rot_Q0':
            {
                'bounds': [0.00000001, 1000000.0],
                'priors': ['Uniform', []],
                'spaces': 'Log_Base10',
                'fixed' : None,
                'unit': 'as input',
            },
        'rot_deltaQ':
            {
                'bounds': [0.00000001, 1000000.0],
                'priors': ['Uniform', []],
                'spaces': 'Log_Base10',
                'fixed' : None,
                'unit': 'as input',
            },
        # GP framework parameter
        'Vc':
            {
                'bounds': [-500.0, 500.0],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : None,
                'unit': 'days',
            },
        'Vr':
            {
                'bounds': [-500.0, 500.0],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : None,
                'unit': 'as input',
            },
        'Lc':
            {
                'bounds': [-500.0, 500.0],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : None,
                'unit': 'as input',
            },
        'Bc':
            {
                'bounds': [-500.0, 500.0],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : None,
                'unit': 'as input',
            },
        'Br':
            {
                'bounds': [-500.0, 500.0],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : None,
                'unit': 'as input',
            },
        # Multi-dimensional parameters
        'rot_amp': # the rotational term
            {
                'bounds': [-500.0, 500.0],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : None,
                'unit': 'as input',
            },
        'con_amp': # the convective blueshift suppression term
            {
                'bounds': [-500.0, 500.0],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : None,
                'unit': 'as input',
            },
        'cos_amp': # the rotational term
            {
                'bounds': [-500.0, 500.0],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : None,
                'unit': 'as input',
            },
        'cos_der': # the convective blueshift suppression term
            {
                'bounds': [-500.0, 500.0],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : None,
                'unit': 'as input',
            },
        'cyc_amp': # the rotational term
            {
                'bounds': [-500.0, 500.0],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : None,
                'unit': 'as input',
            },
        'cyc_der': # the convective blueshift suppression term
            {
                'bounds': [-500.0, 500.0],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : None,
                'unit': 'as input',
            },
        # Matern kernel
        'matern32_sigma': # Matern kernel,  sigma parameter (amplitude^2)
            {
                'bounds': [0.000001, 1000000.0],
                'priors': ['Uniform', []],
                'spaces': 'Log_Base10',
                'fixed' : None,
                'unit': 'as input',
            },
        'matern32_rho': # Matern kernel,  rho parameter (metric)
            {
                'bounds': [0.001, 1000.00],
                'priors': ['Uniform', []],
                'spaces': 'Log_Base10',
                'fixed' : None,
                'unit': 'as input',
            },
        'matern32_multigp_sigma': # Matern kernel,  rho parameter (metric)
            {
                'bounds': [-10000.00, 1000.00],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : None,
                'unit': 'as input',
            },
        'matern32_multigp_sigma_deriv': # Matern kernel,  rho parameter (metric)
            {
                'bounds': [-10000.00, 1000.00],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : None,
                'unit': 'as input',
            },
        # sinusoid-like activity
        'sin_P':  # Period
            {
                'bounds': [1.0, 1000.0],
                'priors': ['Uniform', []],
                'spaces': 'Log_Base2',
                'fixed' : None,
                'unit': 'days',
            },
        'sin_K':  # Sinusoid semi-amplitude
            {
                'bounds': [0.001, 2000.0],
                'priors': ['Uniform', []],
                'spaces': 'Log_Base2',
                'fixed' : None,
                'unit': 'as input',
            },
        'sin_f':  # Sinusoid curve phase
            {
                'bounds': [1.0, 1000.0],
                'priors': ['Uniform', []],
                'spaces': 'Log_Base2',
                'fixed' : None,
                'unit': 'degrees',
            },
    }



    for i_k in range(0,10):
        # Following the priors definitions in Barros+2022
        parameters_dictionary['grn_k'+repr(i_k) + '_period'] = {
                'bounds': [1e-08, 100.],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.00,
                'unit': 'days',
        }
        parameters_dictionary['grn_k'+repr(i_k) + '_sigma'] = {
                'bounds': [1e-08, 100.],
                'priors': ['Uniform', []],
                'spaces': 'Log_Base2',
                'fixed' : None,
                'unit': 'as input',
        }


        parameters_dictionary['osc_k'+repr(i_k) + '_period'] = {
                'bounds': [1e-08, 100.],
                'priors': ['Uniform', []],
                'spaces': 'Linear',
                'fixed' : 0.00,
                'unit': 'days',
        }
        parameters_dictionary['osc_k'+repr(i_k) + '_sigma'] = {
                'bounds': [1e-08, 100.],
                'priors': ['Uniform', []],
                'spaces': 'Log_Base2',
                'fixed' : None,
                'unit': 'as input',
        }
        parameters_dictionary['osc_k'+repr(i_k) + '_Q0'] = {
                'bounds': [1.00, 1e3],
                'priors': ['Uniform', []],
                'spaces': 'Log_Base10',
                'fixed' : None,
                'unit': 'as input',
            }

    # legacy parameters for the second implementation of celerite
    ##'Hamp_factor',
    ##'Q0',  # celerite rotation term parameter
    ##'deltaQ', # celerite rotation term parameter
    ##'ln_Q0',  # celerite rotation term parameter
    ##'ln_deltaQ', # celerite rotation term parameter
    ##'mix', # celerite rotation term parameter
    ##'amp', # celerite rotation term parameter
    ### legacy parameters for the first implementation of celerite
    ##'cel_factor',
    ##'cel_a',  # celerite term A
    ##'cel_b',  # celerite term B
    ##'cel_c',  # celerite term C

    """These default boundaries are used when the user does not define them in the yaml file
    default_bounds = {
        'Hamp_factor': [0.01, 10.0],
        'Q0': [0.00001, 10000.000],
        'deltaQ': [0.00001, 10000.000],
        'mix': [0.001, 1000.000],
        'amp': [0.0001, 1000.0],
        #
        'cel_factor': [0.00000001, 1000000.0],
        'cel_B': [0.00000001, 1000000.0],
        'cel_C': [0.00000001, 1000000.0],
    }

    These default priors are used when the user does not define them in the yaml file
    The boundaries are not included in this dictionary, because it is likely that the user will specify his
    preferred boundaries without changing the priors

    default_priors = {
        'f': ['Uniform', []],
        #
        'Hamp_factor': ['Uniform', []],
        'Q0': ['Uniform', []],
        'deltaQ': ['Uniform', []],
        'mix': ['Uniform', []],
        'amp': ['Uniform', []],
        #
        'cel_factor': ['Uniform', []],
        'cel_B': ['Uniform', []],
        'cel_C': ['Uniform', []],
    }

    default_spaces = {
        'Hamp_factor': 'Linear',
        'Q0': 'Log_Base10',
        'deltaQ': 'Log_Base10',
        'mix': 'Log_Base10',
        'amp': 'Log_Base2',
        #
        'cel_factor': 'Linear',
        'cel_B': 'Logarithmic',
        'cel_C': 'Logarithmic',
    }
    """

    default_fixed = {}

    recenter_pams = {'f'}

    def initialize_model(self, mc, **kwargs):

        self.use_stellar_rotation_period = False
        for keyword in keywords_stellar_rotation:
            self.use_stellar_rotation_period = kwargs.get(keyword, self.use_stellar_rotation_period)

        self.use_stellar_activity_decay = False
        for keyword in keywords_stellar_activity_decay:
            self.use_stellar_activity_decay = kwargs.get(keyword, self.use_stellar_activity_decay)
