import os
import sys


if 'celerite' not in sys.modules:

    try:
        import celerite
        import autograd.numpy as np
    except:
        from dummy import celerite
        import numpy as np
        print('WARNING! Imported dummy celerite, models relying on this package will not work')

    try:
        import PyPolyChord as PyPolyChord
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


def giveback_priors(kind, pams, val):
    if kind == 'Gaussian':
        return -(val - pams[0]) ** 2 / (2 * pams[1] ** 2)
    if kind == 'Uniform':
        return 0.0


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

