import numpy as np
import kepler_exo as kp
import yaml
import george
import spotpy
import ttvfast
import os
import sys
import gc

if os.path.isdir('/Users/malavolta/Astro/CODE/'):
    sys.path.insert(0, '/Users/malavolta/Astro/CODE/others/PolyChord/')
else:
    sys.path.insert(0, '/home/malavolta/CODE/others/PolyChord/')
import PyPolyChord.PyPolyChord as PolyChord

if os.path.isdir('/Users/malavolta/Astro/CODE/'):
    sys.path.insert(0, '/Users/malavolta/Astro/CODE/trades/pytrades/')
else:
    sys.path.insert(0, '/home/malavolta/CODE/trades/pytrades/')
#from pytrades_lib import pytrades
from dummy import pytrades
from dummy import PolyChord
import constants


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
