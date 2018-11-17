import os
import sys
from scipy import stats
from scipy.interpolate import interp1d

if 'celerite' not in sys.modules:

    try:
        from pyde.de import DiffEvol
    except:
        from dummy import DiffEvol
        print('WARNING! Imported dummy pyde, nothing will work')

    try:
        import celerite
        import autograd.numpy as np
    except:
        import numpy as np
        from dummy import Celerite_QuasiPeriodicActivity
        print('WARNING! Imported dummy celerite, models relying on this package will not work')
    
    try:
        import PyPolyChord
        from PyPolyChord.settings import PolyChordSettings
    except:
        from dummy import PyPolyChord
        print('WARNING! Imported dummy PyPolyChord, models relying on this package will not work')

    try:
        if os.path.isdir('/Users/malavolta/Astro/CODE/'):
            sys.path.insert(0, '/Users/malavolta/Astro/CODE/trades/pytrades/')
        else:
            sys.path.insert(0, '/home/malavolta/CODE/trades/pytrades/')
        from pytrades_lib import pytrades
    except:
        from dummy import pytrades
        print('WARNING! Imported dummy TRADES, models relying on this package will not work')

    try:
        import ttvfast
    except:
        from dummy import ttvfast
        print('WARNING! Imported dummy TTVFAST, models relying on this package will not work')

    try:
        import george
    except:
        from dummy import george
        print('WARNING! Imported dummy george, models relying on this package will not work')

import kepler_exo
import yaml
import constants
import gc

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
    return np.sqrt(np.square(ecoso, dtype=np.double) + np.square(esino, dtype=np.double))


def get_2var_e(var, fix, i):
    if len(np.shape(var)) == 1:
        ecoso = var[i[0]]
        esino = var[i[1]]
    else:
        ecoso = var[:, i[0]]
        esino = var[:, i[1]]
    return np.square(ecoso, dtype=np.double) + np.square(esino, dtype=np.double)


def get_2var_o(var, fix, i):
    if len(np.shape(var)) == 1:
        ecoso = var[i[0]]
        esino = var[i[1]]
    else:
        ecoso = var[:, i[0]]
        esino = var[:, i[1]]
    return np.arctan2(esino, ecoso, dtype=np.double)


def get_2darray_from_val(val):
    out = np.zeros(2, dtype=np.double)
    try:
        if len(np.shape(val)) == 1:
            out[0] = val[0]
        else:
            out[:] = val[0:2]
    except:
        out[0] = val
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
Special subroutine to transform MultiNest/PolyChord priors, i.e., trasnform the datacube from [0:1] to physical 
values while taking into account the priors
"""

def nested_sampling_prior_transformation(kind, bounds, pams):

    #x = np.linspace(0.000000, 1.000000, num=1001, endpoint=True)
    x_var = np.linspace(0.000000, 1.000000, num=10001, endpoint=True, dtype=np.double)*(bounds[1]-bounds[0]) + bounds[0]
    area = np.zeros(len(x_var), dtype=np.double)

    for x_num, x_val in enumerate(x_var):
        area[x_num:] += np.exp(giveback_priors(kind, bounds, pams, x_val)) * (1. / 10000.) + 0.000000000001
    area[0] = 0
    area /= area[-1]

    return interp1d(area, x_var, kind='cubic')


def nested_sampling_prior_prepare(kind, bounds, pams):
    """
    This subroutine computes the coefficient of the spline interpolation of the inverse cumulative function
    In some special cases, ruterns the parameters required by the intrinsic function, e.g. scipi.stats.norm.icf
    according to their implementation in nested_sampling_prior_compute()

    :param kind:
    :param bounds:
    :param pams:
    :return:
    """
    
    if kind == 'Uniform' or kind=='Gaussian':
        return pams


def compute_value_sigma(samples):
    if np.size(np.shape(samples)) == 1:
        sample_med = np.zeros(3)
        sample_tmp = np.percentile(samples, [15.865, 50, 84.135], axis=0)
        sample_med[0] = sample_tmp[1]
        sample_med[1] = sample_tmp[2] - sample_tmp[1]
        sample_med[2] = sample_tmp[0] - sample_tmp[1]

    elif np.size(np.shape(samples)) == 2:
        sample_med = np.asarray(map(lambda v: (v[1], v[2] - v[1], v[1] - v[0]),
                                   zip(*np.percentile(samples, [15.865, 50, 84.135], axis=0))))

    else:
        print 'ERROR!!! '
        return None
    return sample_med


def pick_MAP_parameters(samples, lnprob):

    indmax = np.argmax(lnprob)

    if np.size(np.shape(samples)) == 1:
        return samples[indmax], lnprob[indmax]
    elif np.size(np.shape(samples)) == 2:
        return samples[indmax, :], lnprob[indmax]
    else:
        print 'ERROR!!! '
        return None
