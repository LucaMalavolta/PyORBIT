from classes.model_container import ModelContainer
from classes.input_parser import yaml_parser
import classes.kepler_exo as kp
import numpy as np
import h5py
import cPickle as pickle
import scipy.optimize
import csv
import os
import argparse
import matplotlib as mpl
from matplotlib.ticker import FormatStrFormatter
import sys
mpl.use('Agg')
from matplotlib import pyplot as plt
import corner
sys.path.append('/Users/malavolta/Astro/CODE/trades/pytrades')
import constants

#Plot improvements
plt.rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
plt.rc('text', usetex=True)


def GelmanRubin(chains_T):
    # Courtesy of Luca "Sbuffo" Borsato
    n, M = np.shape(chains_T)

    theta_m = [np.mean(chains_T[:,i_m]) for i_m in range(0, M)]
    theta = np.mean(theta_m)

    d_theta2 = (theta_m - theta)**2
    B_n = np.sum(d_theta2) / (M-1)

    arg_W = [np.sum((chains_T[:,i_m] - theta_m[i_m])**2) / (n-1) for i_m in range(0, M)]
    W = np.mean(arg_W)

    n_frac = (n-1)/n
    var_plus = n_frac*W + B_n
    Var = var_plus + (B_n/M)

    Rc = np.sqrt(Var / W)
    return Rc


def get_mass(M_star2, M_star1, Period, K1, e0):
    # M_star1, M_star2 in solar masses
    # P in days -> Period is converted in seconds in the routine
    # inclination assumed to be 90 degrees
    # Gravitational constant is given in m^3 kg^-1 s^-2
    # output in m/s
    output = K1 - (2. * np.pi * G_grav * M_sun / 86400.0) ** (1.0 / 3.0) * (1.000 / np.sqrt(1.0 - e0 ** 2.0)) * (
                                                                                                                    Period) ** (
                                                                                                                    -1.0 / 3.0) * (
                      M_star2 * (M_star1 + M_star2) ** (-2.0 / 3.0))
    return output

G_grav = constants.Gsi # Gravitational Constants in SI system [m^3/kg/s^2]
G_ttvfast = constants.Giau  # G [AU^3/Msun/d^2]
M_SJratio = constants.Msjup
M_SEratio = constants.Msear
M_JEratio = constants.Mjear

R_SJratio = constants.Rsjup
R_JEratio = constants.Rjear
R_SEratio = constants.Rsjup * constants.Rjear

M_sun = constants.Msun
M_jup = constants.Mjup

Mu_sun = constants.Gsi * constants.Msun
seconds_in_day = constants.d2s
AU_km = constants.AU
AUday2ms = AU_km / seconds_in_day * 1000.0


parser = argparse.ArgumentParser(prog='PyORBIT_V4_get_chains.py', description='Extract results from output MCMC')
parser.add_argument('config_file', type=str, nargs=1, help='config file')

args = parser.parse_args()

sampler = args.sample[0]
file_conf = args.config_file[0]

mc = ModelContainer()
yaml_parser(file_conf, mc)

mc.initialize_model()

if bool(mc.pcv.dynamical):
    mc.dynamical_model.prepare(mc, mc.pcv)

M_star1 = mc.star_mass[0]
M_star1_err = mc.star_mass[1]


dir_input = './' + mc.planet_name + '/emcee/'
dir_output = './' + mc.planet_name + '/emcee_plot/'
os.system('mkdir -p ' + dir_output)

mc.variable_list = pickle.load(open(dir_input + 'vlist.pick', 'rb'))
mc.scv.use_offset = pickle.load(open(dir_input + 'scv_offset.pick', 'rb'))

print mc.variable_list

h5f = h5py.File(dir_input + mc.planet_name + '.hdf5', "r")

h5f_data = h5f['/data']
h5f_emcee = h5f['/emcee']

for item in h5f_emcee.attrs.keys():
    print item + ":", h5f_emcee.attrs[item]

    if item == 'nwalkers': mc.emcee_parameters['nwalkers']  = h5f_emcee.attrs[item]
    if item == 'ndim': mc.ndim = h5f_emcee.attrs[item]

mc.bounds = h5f['/emcee/bound']
chain = h5f['/emcee/chain']
lnprobability = h5f['/emcee/lnprobability']
acceptance_fraction = h5f['/emcee/acceptance_fraction']
#acor = h5f['/emcee/acor']

print mc.bounds[:]

print
print '*************************************************************'
print
print 'Acceptance Fraction for all walkers:'
print acceptance_fraction[:]
print
print '*************************************************************'
print

if args.nburn > 0:
    mc.emcee_parameters['nburn'] = args.nburn

if mc.emcee_parameters['nsave'] > 0:
    mc.emcee_parameters['nsteps'] = h5f_emcee.attrs['nsample']
    if mc.emcee_parameters['nburn'] > mc.emcee_parameters['nsteps']:
        mc.emcee_parameters['nburn'] = mc.emcee_parameters['nsteps'] / 4

ntotal = np.int(mc.emcee_parameters['nsteps'] / mc.emcee_parameters['thin'])
nburnin = np.int(mc.emcee_parameters['nburn'] / mc.emcee_parameters['thin'])

lnprb_T = lnprobability[:][:].T

chain_T = np.ndarray([ntotal, mc.emcee_parameters['nwalkers'], mc.ndim], dtype=np.double)
for ii in xrange(0, mc.ndim):
    chain_T[:, :, ii] = chain[:, :, ii].T

'''chain_T contains the chains coming out from emcee, after being transposed 
the individual chain for walker m of parameter n is given by  chain_T[:,]

'''

chain_burnt = chain_T[nburnin:, :, :]
s = chain_burnt.shape
lnprob_burnt = lnprb_T[nburnin:, :,]
flatchain = chain_burnt.reshape(s[1] * s[0], s[2])
flatlnprob = lnprob_burnt.reshape(s[1] * s[0])
nsample = s[1] * s[0]
#sel_flatchain = flatchain[:, 0] < 1.

chain_med = np.asarray(map(lambda v: (v[1], v[2] - v[1], v[1] - v[0]),
                           zip(*np.percentile(flatchain[:, :], [15.865, 50, 84.135], axis=0))))
lnprob_med = np.percentile(flatlnprob, [15.865, 50, 84.135], axis=0)
lnprob_med[1:] = np.abs(lnprob_med[1:]-lnprob_med[0])

mc.results_resumen(chain_med[:, 0])

print
print '*************************************************************'
print

lnprob_median = np.median(flatlnprob)
fig = plt.figure(figsize=(12, 12))
plt.xlabel('$\ln \mathcal{L}$')
plt.plot(lnprb_T[:, :], '-', alpha=0.5)
plt.axhline(lnprob_median)
plt.axvline(nburnin, c='r')
plt.savefig(dir_output + 'LNprob_chain.png', bbox_inches='tight', dpi=300)
plt.close(fig)
print 'LNprob median = ', lnprob_median

print
print 'Reference Time Tref: ', mc.Tref
print
print '*************************************************************'
print

# Plotting out the chains
if args.c != 'False':
    for ii in xrange(0, mc.ndim):
        print mc.pam_names[ii], chain_med[ii, 0], ' +\sigma ', chain_med[ii, 1], ' -\sigma ', chain_med[ii, 2]
        file_name = dir_output + 'chain_' + repr(ii) + '_' + mc.pam_names[ii] + '.png'
        fig = plt.figure(figsize=(12, 12))
        plt.plot(chain_T[:, :, ii], '-', alpha=0.5)
        plt.axvline(nburnin, c='r')
        plt.axhline(chain_med[ii, 0], c='k')
        plt.savefig(file_name, bbox_inches='tight', dpi=300)
        plt.close(fig)

print
print '*************************************************************'
print

if args.t != 'False':
    for nd in xrange(0, mc.ndim):  # (0,ndim):
        out_absc = np.arange(0, nburnin, 1)
        out_lines = np.zeros(nburnin)
        for ii in xrange(20, nburnin):
            out_lines[ii] = GelmanRubin(chain_T[:ii, :, nd])

        fig = plt.figure(figsize=(12, 12))
        #plt.ylim(0.95, 2.3)
        plt.plot(out_absc[20:], out_lines[20:], '-', color='k')
        plt.axhline(1.01)
        plt.savefig(dir_output + 'GRtrace_pam_' + repr(nd) + '.png', bbox_inches='tight', dpi=300)
        plt.close(fig)

    print
    print '*************************************************************'
    print


""" Creation of the directory for the plots"""
plot_dir = dir_output + '/files_plot/'

if not os.path.exists(plot_dir):
    os.makedirs(plot_dir)

boundaries = np.asarray([mc.Tref, mc.Tref])
plot_dir = dir_output + '/files_plot/'

if not os.path.exists(plot_dir):
    os.makedirs(plot_dir)
for dataset in mc.dataset_dict.itervalues():
    if dataset.kind == 'RV':
        boundaries[0] = min(boundaries[0], dataset.x[0])
        boundaries[1] = max(boundaries[1], dataset.x[-1])

boundaries += np.asarray([-1.0, 1.0]) * (boundaries[1] - boundaries[0]) / 10.

x_range_step = max(0.01, (boundaries[1] - boundaries[0]) / 100000)
x_range = np.arange(boundaries[0], boundaries[1], x_range_step)
x_phase = np.arange(-0.50, 1.50, 0.005, dtype=np.double)

model_dsys, model_plan, model_orbs, model_actv, model_curv = mc.rv_make_model(chain_med[:, 0], x_range, x_phase)

