import os
import sys

sys.path.insert(0, '/Users/malavolta/Astro/CODE/pyorbit/')
import pyorbit
import pyorbit.classes.kepler_exo as kp

# creating




#input_dataset{'RV'} =

config_in = {
    'output': 'test0_1planet_dataset_V5',
    'inputs': {
        'RV': {
            'models': ['radial_velocities'],
            'kind': 'RV',
            'file': 'test_1planet_RV.dat'
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
        'star_radius': [1.0, 0.01]
    },
    'solver': {
        'pyde': {
            'npop_mult': 4,
            'ngen': 8000
        },
        'emcee': {
            'npop_mult': 4,
            'nburn': 20000,
            'nsteps': 50000,
            'thin': 100
        },
        'recenter_bounds': True
    }
}

pyorbit.pyorbit_emcee(config_in)