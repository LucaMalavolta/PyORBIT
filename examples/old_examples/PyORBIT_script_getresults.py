import matplotlib
matplotlib.use('WXAgg')
import matplotlib.pyplot as plt


import os
import sys
import numpy as np

sys.path.insert(0, '/home/malavolta/CODE/pyorbit/')
sys.path.insert(0, '/Users/malavolta/Astro/CODE/pyorbit/')
import pyorbit
import pyorbit.classes.kepler_exo as kp


"""This is the config file as it would appear after being read by the terminal interface
Two things to be noticed:
    1) the keyword 'output' must be explicited (in the terminal version, it's created starting from the file name)
    2) the 'file' keyword has been commented out since the dataset is provided internally to the script
Notice that the 'file ' keyword has been commented, because we are providing
"""

config_in = {
    'output': 'test0_1planet_dataset_V5',
    'inputs': {
        'RV': {
            'models': ['radial_velocities'],
            'kind': 'RV',
            #'file': 'test_1planet_RV.dat'
        }
    },
    'models': {
        'radial_velocities': {
            'planets': ['b']
        }
    },
    'common': {
        'planets': {
            'b': {
                'boundaries': {
                    'P': [2.0, 100.0],
                    'K': [0.01, 100.0],
                    'e': [0.0, 1.0]
                },
                'fixed': {
                    'i': [90.0, 0.001]
                },
                'orbit': 'keplerian'
            }
        }
    },
    'parameters': {
        'star_mass': [1.0, 0.01],
        'star_radius': [1.0, 0.01],
    },
    'solver': {
        'pyde': {
            'npop_mult': 4,
            'ngen': 2000
        },
        'emcee': {
            'npop_mult': 4,
            'nburn': 3000,
            'nsteps': 10000,
            'thin': 100
        },
        'recenter_bounds': True
    }
}



"""
Set to true the plots you want
"""
plot_dictionary = {
    'chains': False,
    'plot': False,
    'traces': False,
    'full_correlation': True,
    'common_corner': True,
    'dataset_corner': False,
    'model_files': True,
    'MAP_model_files': False
}
sampler = 'emcee'

pyorbit.pyorbit_getresults(config_in, sampler, plot_dictionary)
