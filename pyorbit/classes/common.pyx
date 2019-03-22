from __future__ import print_function
import os
import sys
from scipy import stats
from scipy.interpolate import splrep, splev
import pyorbit.classes.constants as constants
import yaml
import gc

if 'celerite' not in sys.modules:

    try:
        import dynesty
    except:
        from pyorbit.classes.dummy import dynesty
        print('WARNING! Imported pyorbit.classes.dummy dynesty, models relying on this package will not work')

    try:
        from pyde.de import DiffEvol
    except:
        from pyorbit.classes.dummy import DiffEvol
        print('WARNING! Imported pyorbit.classes.dummy pyde, nothing will work')

    try:
        import celerite
        import autograd.numpy as np
    except:
        import numpy as np
        from pyorbit.classes.dummy import Celerite_QuasiPeriodicActivity
        print('WARNING! Imported pyorbit.classes.dummy celerite, models relying on this package will not work')
    
    try:
        import PyPolyChord
        from PyPolyChord.settings import PolyChordSettings
    except:
        from pyorbit.classes.dummy import PyPolyChord
        print('WARNING! Imported pyorbit.classes.dummy PyPolyChord, models relying on this package will not work')

    try:
        if os.path.isdir('/Users/malavolta/Astro/CODE/'):
            sys.path.insert(0, '/Users/malavolta/Astro/CODE/trades/pytrades/')
        else:
            sys.path.insert(0, '/home/malavolta/CODE/trades/pytrades/')
        from pytrades_lib import pytrades
    except:
        from pyorbit.classes.dummy import pytrades
        print('WARNING! Imported pyorbit.classes.dummy TRADES, models relying on this package will not work')

    try:
        import ttvfast
    except:
        from pyorbit.classes.dummy import ttvfast
        print('WARNING! Imported pyorbit.classes.dummy TTVFAST, models relying on this package will not work')

    try:
        import george
    except:
        from pyorbit.classes.dummy import george
        print('WARNING! Imported pyorbit.classes.dummy george, models relying on this package will not work')

    try:
        import batman
    except:
        from pyorbit.classes.dummy import batman
        print('WARNING! Imported pyorbit.classes.dummy batman, models relying on this package will not work')

def get_var_log(var, fix, i):
    if len(np.shape(var)) == 1:
        return np.log2(var[i], dtype=np.double)
    else:
        return np.log2(var[:, i], dtype=np.double)


def get_var_exp(var, fix, i):
    if len(np.shape(var)) == 1:
        return np.exp2(var[i], dtype=np.double)
    else:
        return np.exp2(var[:, i], dtype=np.double)


def get_var_val(var, fix, i):
    if len(np.shape(var)) == 1:
        return var[i]
    else:
        return var[:, i]


def get_fix_val(var, fix, i):
    if len(np.shape(fix)) == 1:
        return fix[i]
    else:
        return fix[:, i]


def get_2var_sre(var, fix, i):
    if len(np.shape(var)) == 1:
        ecoso = var[i[0]]
        esino = var[i[1]]
    else:
        ecoso = var[:, i[0]]
        esino = var[:, i[1]]
    return np.square(ecoso, dtype=np.double) + np.square(esino, dtype=np.double)


def get_2var_e(var, fix, i):
    if len(np.shape(var)) == 1:
        ecoso = var[i[0]]
        esino = var[i[1]]
    else:
        ecoso = var[:, i[0]]
        esino = var[:, i[1]]
    return np.sqrt(np.square(ecoso, dtype=np.double) + np.square(esino, dtype=np.double))


def get_2var_o(var, fix, i):
    if len(np.shape(var)) == 1:
        ecoso = var[i[0]]
        esino = var[i[1]]
    else:
        ecoso = var[:, i[0]]
        esino = var[:, i[1]]
    return np.arctan2(esino, ecoso, dtype=np.double)


# Kipping transformation for quadratic limb darkening coefficients
def get_2var_c1(var, fix, i):
    if len(np.shape(var)) == 1:
        q1 = var[i[0]]
        q2 = var[i[1]]
    else:
        q1 = var[:, i[0]]
        q2 = var[:, i[1]]
    return 2.0*np.sqrt(q1, dtype=np.double) * q2


def get_2var_c2(var, fix, i):
    if len(np.shape(var)) == 1:
        q1 = var[i[0]]
        q2 = var[i[1]]
    else:
        q1 = var[:, i[0]]
        q2 = var[:, i[1]]
    return np.sqrt(q1, dtype=np.double) * (1.0 - 2.0*q2)


def get_2darray_from_val(val):
    out = np.zeros(2, dtype=np.double)
    try:
        if np.shape(val) == 1:
            out[0] = val[0]
        else:
            out[:] = val[0:2]
    except:
        try:
            out[0] = val
        except:
            out = val
    return out


def giveback_priors(kind, bounds, pams, val):

    """ The code is supposed to give -np.inf log-likelihood when the parameters are outside the boundaries,
    so this case is not emcopassed in the definition of the priors """

    if kind == 'None':
        return 0.00

    if kind == 'Gaussian':
        return -(val - pams[0]) ** 2 / (2 * pams[1] ** 2)

    if kind == 'Uniform':
        return np.log(1./(bounds[1]-bounds[0]))

    if kind == 'Jeffreys' or kind == 'TruncatedJeffreys':
        return np.log(1./(val*np.log(bounds[1]/bounds[0])))

    if kind == 'ModifiedJeffreys' or kind == 'TruncatedModifiedJeffreys':
        """ Used for the semi-amplitude of the RV curve
            bounds[1] = Kmax (suggested 999 m/s)
            pams[0] = K_0 (suggested 1 m/s)
        """
        return np.log(1./(pams[0]*(1. + val/pams[0])) * 1./np.log(1.+bounds[1]/pams[0]))

    if kind == "TruncatedRayleigh":
        """ bounds[1] = e_max
            pams[0] = sigma_2 (suggested 0.2)
        """
        return np.log((val/pams[0]**2 * np.exp(-0.5 *(val/pams[0])**2))/(1. - np.exp(-0.5*(bounds[1]/pams[0])**2)))

    if kind == "WhiteNoisePrior":
        """ bounds[1] = noise_max (99 m/s)
            pams[0] = noise_0 (suggested 1 m/s)
        """
        return np.log(1./(pams[0]*(1.0 + val/pams[0])) * 1.0/np.log(1.0+bounds[1]/pams[0]))

    if kind == "BetaDistribution":
        return np.log(stats.beta.pdf((val-bounds[0])/(bounds[1]-bounds[0]), pams[0], pams[1]))


"""
DEPRECATED
Special subroutine to transform MultiNest/PolyChord priors, i.e., trasnform the datacube from [0:1] to physical 
values while taking into account the priors

def nested_sampling_prior_transformation(kind, bounds, pams):

    #x = np.linspace(0.000000, 1.000000, num=1001, endpoint=True)
    x_var = np.linspace(0.000000, 1.000000, num=10001, endpoint=True, dtype=np.double)*(bounds[1]-bounds[0]) + bounds[0]
    area = np.zeros(len(x_var), dtype=np.double)

    for x_num, x_val in enumerate(x_var):
        area[x_num:] += np.exp(giveback_priors(kind, bounds, pams, x_val)) * (1. / 10000.) + 0.000000000001
    area[0] = 0
    area /= area[-1]

    return interp1d(area, x_var, kind='cubic')
"""


def nested_sampling_prior_prepare(kind, bounds, pams, space):
    """
    This subroutine computes the coefficient of the spline interpolation of the inverse cumulative function
    In some special cases, ruterns the parameters required by the intrinsic function, e.g. scipi.stats.norm.icf
    according to their implementation in nested_sampling_prior_compute()

    :param kind: type of prior
    :param bounds: list/array with lower and upper limits for parameter exploration
    :param pams: parameters relative to the prior probability function
    :return:
    """
    
    if kind == 'Uniform':
        return bounds

    if kind == 'Gaussian':
        return pams

    if kind == 'beta':
        return pams

    """ All the following priors are defined only if the variable is sampled in the Natural space"""
    if space is not 'Linear':
        print()
        print(' *** ERROR in the YAML file ***')
        print(' You are using a prior that is not supported in a not-Linear sampling space')
        print(' add this keyword in the YAML file for each parameter not sampled in the ')
        print(' Linear space and with a prior other than Uniform or Gaussian')
        print('   spaces: ')
        print('       pam: Linear')
        print()
        quit()

    x_var = np.linspace(0.000000, 1.000000, num=10001, endpoint=True, dtype=np.double)*(bounds[1]-bounds[0]) + bounds[0]
    area = np.zeros(len(x_var), dtype=np.double)

    for x_num, x_val in enumerate(x_var):
        area[x_num:] += np.exp(giveback_priors(kind, bounds, pams, x_val)) * (1. / 10000.) + 0.000000000001
    area[0] = 0
    area /= area[-1]

    return splrep(area, x_var)


def nested_sampling_prior_compute(val, kind, coeff, space):
    """
    In same cases ()

    :param val:
    :param kind:
    :param coeff:
    :param space:
    :return:
    """

    if kind == 'Uniform':
        return val * (coeff[1] - coeff[0]) + coeff[0]

    if kind == 'Gaussian':
        if space == 'Logarithmic':
            return np.log2(stats.norm.isf(val, coeff[0], coeff[1]))
        else:
            return stats.norm.isf(val, coeff[0], coeff[1])

    elif kind == 'beta':
        if space == 'Logarithmic':
            return np.log2(stats.norm.isf(stats.beta.isf(val, coeff[0], coeff[1])))
        else:
            return stats.beta.isf(val, coeff[0], coeff[1])

    return splev(val, coeff)


def compute_value_sigma(samples):
    if np.size(np.shape(samples)) == 1:
        sample_med = np.zeros(3)
        sample_tmp = np.percentile(samples, [15.865, 50, 84.135], axis=0)
        sample_med[0] = sample_tmp[1]
        sample_med[1] = sample_tmp[2] - sample_tmp[1]
        sample_med[2] = sample_tmp[0] - sample_tmp[1]

    elif np.size(np.shape(samples)) == 2:
        sample_med = np.asarray(list(map(lambda v: (v[1], v[2] - v[1], v[1] - v[0]),
                                   zip(*np.percentile(samples, [15.865, 50, 84.135], axis=0)))))

    else:
        print('ERROR!!! ')
        return None
    return sample_med


def pick_MAP_parameters(samples, lnprob):

    indmax = np.argmax(lnprob)
    if np.size(np.shape(samples)) == 1:
        return samples[indmax], lnprob[indmax]
    elif np.size(np.shape(samples)) == 2:
        return samples[indmax, :], lnprob[indmax]
    else:
        print('ERROR!!! ')
        return None


def convert_rho_to_a(P, rho):

    return np.power(constants.Gsi * (constants.d2s * constants.d2s) * (P**2)
                    * rho * constants.rho_Sun / (3. * np.pi), 1./3.)


def convert_b_to_i(b,e,o,a):

    rho_e = (1. - e ** 2) / (1. + e * np.sin(o))
    arccos_argument = b / a / rho_e
    if np.size(arccos_argument)<=1:

        if arccos_argument > 1.: arccos_argument=1
        if arccos_argument < -1.: arccos_argument=-1
    else:
        arccos_argument = [1. if b > 1. else b for b in arccos_argument]
        arccos_argument = [-1. if b < -1. else b for b in arccos_argument]

    return np.arccos(arccos_argument)*180./np.pi


