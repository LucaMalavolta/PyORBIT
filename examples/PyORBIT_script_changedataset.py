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

# creating

""" Fake input dataset is created"""
x = np.random.normal(np.arange(6000, 6100, 1, dtype=np.double), 0.2)
"Let's avoid 1-day alias"
n_x = np.size(x)

Tref = 6049.48253491
x0 = x - Tref

P = 13.4237346
K = 23.47672
phase = 0.34658203
e = 0.13
omega = 0.673434
offset = 45605

""" These is the syntethic dataset without error """
y_pla = kp.kepler_RV_T0P(x0, phase, P, K, e, omega) + offset

""" Array with the associated error"""
rv_err = np.ones(n_x) * 2.0

""" adding white noise to the data"""
mod_pl = np.random.normal(y_pla, rv_err)

"""When provided to the pyorbit_emcee() subroutine instead of being read from a file, the input_dataset must be a
dictionary where the keyword corresponds to the label of the dataset. for each keyword, a n*6 numpy array must be
included, where n is the number of observations.
If the keyword for a dataset is present, it will have priority on the dataset file unless the keyword is empty

For example: input_dataset{'RV'} = np.zeros([n,6])
required:
input_dataset{'RV'}[:,0] = time of observation
input_dataset{'RV'}[:,1] = value of observation
input_dataset{'RV'}[:,2] = associated error
input_dataset{'RV'}[:,4] = jitter flag (stars from 0, default is -1, no jitter)
input_dataset{'RV'}[:,5] = offset flag (stars from 0, default is -1, no offset)
input_dataset{'RV'}[:,6] = linear trend flag (stars from 0, default is -1, no linear)

The structure is slightly different for Tcent data:
required:
input_dataset{'RV'}[:,0] = number of transit (e.g. if some transit is missing)
input_dataset{'RV'}[:,1] = transit time
input_dataset{'RV'}[:,2] = associated error
input_dataset{'RV'}[:,4] = jitter flag (stars from 0, default is -1, no jitter)
input_dataset{'RV'}[:,5] = offset flag (stars from 0, default is -1, no offset)
input_dataset{'RV'}[:,6] = linear trend flag (stars from 0, default is -1, no linear)
"""

input_dataset = {'RV': np.zeros([n_x, 6])}
input_dataset['RV'][:, 0] = x
input_dataset['RV'][:, 1] = mod_pl
input_dataset['RV'][:, 2] = rv_err
input_dataset['RV'][:, 5] = -1

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
        'Tref:': Tref
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

""" PyORBIT emcee is run using the callable interface
Note the two optional keywords:
    'input_datasets' : can be omitted if all the datasets are read from files
    'return_output' : it will give back sampler_population and sampler_lnprobability, either form compatation or from
    saved files
"""

mc, population, prob = pyorbit.pyorbit_emcee(config_in, input_datasets=input_dataset, return_output=True)

""" Retrieving the emcee parameters used to perform the calculation """
nburnin = mc.emcee_parameters['nburn']
nthin = mc.emcee_parameters['thin']

""" burn-inphase is revomed from each chain (remember: there are _nwalkers_ chains for each parameter)
then the chains are flattened - it means that all the chains for a given parameter are mixed together
Same work for the log-probability array -

"""

flat_chain = pyorbit.classes.io_subroutines.emcee_flatchain(population, nburnin, nthin)
flat_lnprob = pyorbit.classes.io_subroutines.emcee_flatlnprob(prob, nburnin, nthin)

""" The compute_value_sigma will reutrn ndim*3 array with
med[0]: 50th percentile (median)
med[1]: 84.135th percentile - median (+1 sigma)
med[2]: 15.865th percentile - median (-1 sigma)
"""
chain_med = pyorbit.classes.common.compute_value_sigma(flat_chain)
lnprob_med = pyorbit.classes.common.compute_value_sigma(flat_lnprob)
print ' LN probability: %12f   %12f %12f (15-84 p) ' % (lnprob_med[0], lnprob_med[2], lnprob_med[1])

""" Picking the moidel with the maximum loglikelihood"""
theta_maxln = flat_chain[np.argmax(flat_lnprob), :]

"""Let's compute more points from the model"""
x_additional = np.random.normal(np.arange(6120, 6131, 1, dtype=np.double), 0.2)
n_x_additional = np.size(x_additional)
y_pla_additional = kp.kepler_RV_T0P(x_additional-mc.Tref, phase, P, K, e, omega) + offset
rv_err_additional = np.ones(n_x_additional) * 2.0
mod_pl_additional = np.random.normal(y_pla_additional, rv_err_additional)

input_dataset_additional = {'RV': np.zeros([n_x+n_x_additional, 6])}

input_dataset_additional['RV'][:n_x, :] = input_dataset['RV'][:,:]
input_dataset_additional['RV'][n_x:, 0] = x_additional
input_dataset_additional['RV'][n_x:, 1] = mod_pl_additional
input_dataset_additional['RV'][n_x:, 2] = rv_err_additional
input_dataset_additional['RV'][n_x:, 5] = -1

""" Load the new dataset inside the Model Container object
INTERFACE MUST (and it will) IMPROVED!
"""
print 'Outcome BEFORE dataset change '
print 'LN_likelihood of the solution with the median of the posteriors', mc(chain_med[:, 0])
print 'LN_likelihood of the solution with the maximum ln_likelihood', mc(theta_maxln)
print

for key in input_dataset_additional.iterkeys():
    mc.dataset_dict[key].define_dataset_base(input_dataset_additional[key], update=True)

mc.initialize_logchi2()
"""print new ln-likelihood"""
print 'Outcome AFTER dataset change '
print 'LN_likelihood of the solution with the median of the posteriors', mc(chain_med[:, 0])
print 'LN_likelihood of the solution with the maximum ln_likelihood', mc(theta_maxln)
print

""" Plot stuff, to be moved inside a subroutine """

bjd_plot = {
    'full': {
        'start':None, 'end': None, 'range': None
    }
}

kinds = {}
for dataset_name, dataset in mc.dataset_dict.iteritems():
    if dataset.kind in kinds.keys():
        kinds[dataset.kind].extend([dataset_name])
    else:
        kinds[dataset.kind] = [dataset_name]

    bjd_plot[dataset_name] = {
        'start': np.amin(dataset.x0),
        'end': np.amax(dataset.x0),
        'range': np.amax(dataset.x0)-np.amin(dataset.x0),
    }

    if bjd_plot[dataset_name]['range'] < 0.1 : bjd_plot[dataset_name]['range'] = 0.1

    bjd_plot[dataset_name]['start'] -= bjd_plot[dataset_name]['range'] * 0.05
    bjd_plot[dataset_name]['end'] += bjd_plot[dataset_name]['range'] * 0.05
    bjd_plot[dataset_name]['x0_plot'] = \
        np.arange(bjd_plot[dataset_name]['start'], bjd_plot[dataset_name]['end'], 0.1)

    if bjd_plot['full']['range']:
        bjd_plot['full']['start'] = min(bjd_plot['full']['start'], np.amin(dataset.x0))
        bjd_plot['full']['end'] = min(bjd_plot['full']['start'], np.amax(dataset.x0))
        bjd_plot['full']['range'] = bjd_plot['full']['end']-bjd_plot['full']['start']
    else:
        bjd_plot['full']['start'] = np.amin(dataset.x0)
        bjd_plot['full']['end'] = np.amax(dataset.x0)
        bjd_plot['full']['range'] = bjd_plot['full']['end']-bjd_plot['full']['start']

bjd_plot['full']['start'] -= bjd_plot['full']['range']*0.05
bjd_plot['full']['end'] += bjd_plot['full']['range']*0.05
bjd_plot['full']['x0_plot'] = np.arange(bjd_plot['full']['start'],bjd_plot['full']['end'],0.1)

for dataset_name, dataset in mc.dataset_dict.iteritems():
    if dataset.dynamical:
        bjd_plot[dataset_name] = bjd_plot['full']

print 'Plot of the solution using the median of each posterior'
bjd_plot['model_out'], bjd_plot['model_x0'] = mc.get_model(chain_med[:, 0], bjd_plot)

dataset_name = 'RV'
plt.errorbar(mc.dataset_dict[dataset_name].x0,
             mc.dataset_dict[dataset_name].y - bjd_plot['model_out'][dataset_name]['systematics'],
             yerr=mc.dataset_dict[dataset_name]. e,
             fmt='o', zorder=2)
plt.plot(bjd_plot[dataset_name]['x0_plot'], bjd_plot['model_x0'][dataset_name]['complete'], zorder=1)

plt.show()

print 'Plot of the solution with the maximum ln_likelihoodr'
bjd_plot['model_out'], bjd_plot['model_x0'] = mc.get_model(theta_maxln, bjd_plot)

dataset_name = 'RV'
plt.errorbar(mc.dataset_dict[dataset_name].x0,
             mc.dataset_dict[dataset_name].y - bjd_plot['model_out'][dataset_name]['systematics'],
             yerr=mc.dataset_dict[dataset_name].e,
             fmt='o', zorder=2)
plt.plot(bjd_plot[dataset_name]['x0_plot'], bjd_plot['model_x0'][dataset_name]['complete'], zorder=1)

plt.show()

print mc.dataset_dict[dataset_name].x[-2:]
print mc.dataset_dict[dataset_name].x0[-2:]
print mc.dataset_dict[dataset_name].y[-2:]
print bjd_plot['model_out'][dataset_name]['systematics'][-2:]
print mc.dataset_dict[dataset_name].e[-2:]
