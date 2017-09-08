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


parser = argparse.ArgumentParser(prog='PyORBIT_V4_GetResults.py', description='Extract results from output MCMC')
# parser.add_argument('-l', type=str, nargs='+', help='line identificator')
parser.add_argument('sample', type=str, nargs=1, help='sample (emcee or polychord)')
parser.add_argument('config_file', type=str, nargs=1, help='config file')
parser.add_argument('-p', type=str, nargs='?', default='False', help='Create plot files')
parser.add_argument('-mp', type=str, nargs='?', default='False', help='Create MEGA plot - it may be slow!')
parser.add_argument('-v', type=str, nargs='?', default='False', help='Create Veusz ancillary files')
parser.add_argument('-t', type=str, nargs='?', default='False', help='Create Gelman-Rubin traces')
parser.add_argument('-nburn', type=int, nargs='?', default=0, help='Change the emcee burn-in value to the input value')
parser.add_argument('-c', type=str, nargs='?', default='False', help='Create Chains plots')
parser.add_argument('-forecast', type=float, nargs='?', default=None, help='Project RV curve samples into the phase plot')

args = parser.parse_args()

sampler = args.sample[0]
file_conf = args.config_file[0]

sample_keyword = {
    'polychord':['polychord', 'PolyChord', 'polychrod', 'poly'],
    'emcee': ['emcee', 'MCMC', 'Emcee']
}

# file_conf = raw_input()

mc = ModelContainer()
yaml_parser(file_conf, mc)

if sampler in sample_keyword['polychord'] and \
        mc.polychord_parameters['shutdown_jitter']:
    for dataset in mc.dataset_dict.itervalues():
        dataset.shutdown_jitter()

mc.initialize_model()

if bool(mc.models['planets'].dynamical):
    mc.dynamical_model.prepare(mc, mc.models['planets'])

M_star1 = mc.star_mass[0]
M_star1_err = mc.star_mass[1]

if sampler in sample_keyword['emcee']:

    dir_input = './' + mc.planet_name + '/emcee/'
    dir_output = './' + mc.planet_name + '/emcee_plot/'
    os.system('mkdir -p ' + dir_output)

    mc.variable_list = pickle.load(open(dir_input + 'vlist.pick', 'rb'))
    #mc.scv.use_offset = pickle.load(open(dir_input + 'scv_offset.pick', 'rb'))

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

if sampler in sample_keyword['polychord']:

    dir_input = './' + mc.planet_name + '/' + mc.polychord_parameters['base_dir']
    dir_output = './' + mc.planet_name + '/polychord_plot/'
    os.system('mkdir -p ' + dir_output)

    data_in = np.genfromtxt(dir_input+mc.planet_name+'_equal_weights.txt')
    flatlnprob = data_in[:, 1]
    flatchain = data_in[:, 2:]
    nsample = np.size(flatlnprob)

    print
    print '*************************************************************'
    print

    chain_med = np.asarray(map(lambda v: (v[1], v[2] - v[1], v[1] - v[0]),
                               zip(*np.percentile(flatchain[:, :], [15.865, 50, 84.135], axis=0))))
    lnprob_med = np.percentile(flatlnprob, [15.865, 50, 84.135], axis=0)
    lnprob_med[1:] = np.abs(lnprob_med[1:]-lnprob_med[0])
    mc.results_resumen(chain_med[:, 0])

    print
    print '*************************************************************'
    print

    lnprob_median = np.median(flatlnprob)
    print 'LNprob median = ', lnprob_median

    print
    print 'Tref: ', mc.Tref
    print
    print '*************************************************************'
    print

if args.mp != 'False':

    print 'MEGA plot'
    # plotting mega-corner plot
    plt.rc('text', usetex=False)

    megaplot_dat = np.zeros([np.size(flatchain,axis=0), np.size(flatchain, axis=1)+1])
    megaplot_med = np.zeros(np.size(flatchain, axis=1)+1)
    megaplot_dat[:, :-1] = flatchain[:,:]
    megaplot_dat[:, -1] = flatlnprob[:]
    megaplot_med[:-1] = chain_med[:,0]
    megaplot_med[-1] = lnprob_median
    labels = mc.pam_names
    labels.extend(['ln L'])
    fig = corner.corner(megaplot_dat[:, :], labels=labels, truths=megaplot_med)
    fig.savefig(dir_output + "ALL_corners.pdf", bbox_inches='tight', dpi=300)
    plt.close(fig)
    plt.rc('text', usetex=True)

    print
    print '*************************************************************'
    print


#if args.cc != 'False':
#    cc = ChainConsumer()
#    for nd in xrange(0, mc.ndim):  # (0,ndim):
#        cc.add_chain(chain[:, :, nd].flatten(), walkers=mc.nwalkers)
#
#    #print(cc.get_latex_table())
#    print cc.get_summary()
#
#    print cc.diagnostic_gelman_rubin(threshold=0.05)
#    print cc.diagnostic_geweke()
#    print
#    print '*************************************************************'
#    print

x0 = 1. / 150

M_star1_rand = np.random.normal(M_star1, M_star1_err, nsample)

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
# tenpercent = (boundaries[1] - boundaries[0]) / 10.
# boundaries = [boundaries[0] - tenpercent, boundaries[1] + tenpercent]
boundaries += np.asarray([-1.0, 1.0]) * (boundaries[1] - boundaries[0]) / 10.

x_range_step = max(0.01, (boundaries[1] - boundaries[0]) / 100000)
x_range = np.arange(boundaries[0], boundaries[1], x_range_step)
x_phase = np.arange(-0.50, 1.50, 0.005, dtype=np.double)

model_dsys, model_plan, model_orbs, model_actv, model_curv = mc.rv_make_model(chain_med[:, 0], x_range, x_phase)

"""Datasets summary"""

for model in mc.model_list:

    if mc.models[model].model_class == 'planets':

        if sampler in sample_keyword['polychord']:
            ''' Special plot for the polychord case:
                Let's save all the sample together
            '''
            sample_total = {}

        #if args.p != 'False' or args.v != 'False':

        for planet_name in mc.models['planets'].planet_name:

            print 'Planet ', planet_name, ' summary'
            print

            dynamical_flag = (planet_name in mc.models['planets'].dynamical)

            """ Let's feed the conversion function with the average results to get out a list of human variables
                Before doing that, the 1-sigma intervals must be converted to the actual values
            """
            chain_sig = chain_med[:, 0]

            convert_med = mc.models['planets'].convert(planet_name, chain_sig)

            n_orbital = len(mc.models['planets'].var_list[planet_name])
            n_fitted = len(mc.models['planets'].var_list[planet_name]) - len(mc.models['planets'].fix_list[planet_name])

            """ An index is assigned to each variable to keep track of them in
            """
            convert_out = {}
            for n_var, var in enumerate(convert_med):
                convert_out[var] = n_var

            init_var = max(convert_out.values())
            convert_out['Tperi'] = init_var + 1
            convert_out['Tcent'] = init_var + 2
            convert_out['M_kep'] = init_var + 3
            convert_out['a_smj'] = init_var + 4

            if 'curvature' in mc.model_list:
                init_var = max(convert_out.values())
                for n_var, var in enumerate(mc.models[model].list_pams):
                    convert_out[var] = init_var + n_var + 1

            n_total = max(convert_out.values()) + 1 #
            sample_plan = np.zeros([nsample, n_total])
            median_tmp  = np.zeros([n_total])

            """Let's put all the human variables - including those that have been fixed - in sample_plan
               We copy the median and sigma values from the derived distribution"""

            for n_var, var in enumerate(convert_med):
                median_tmp[n_var] = convert_med[var]

            for ii in xrange(0, nsample):
                convert_tmp = mc.models['planets'].convert(planet_name, flatchain[ii, :])
                for var in convert_med:
                    sample_plan[ii, convert_out[var]] = convert_tmp[var]

            # Time of Periastron
            sample_plan[:, convert_out['Tperi']] = mc.Tref + (-sample_plan[:, convert_out['f']] + sample_plan[:, convert_out['o']]) / \
                (2*np.pi) * sample_plan[:, convert_out['P']]

            # Time of transit center
            sample_plan[:, convert_out['Tcent']] = mc.Tref + kp.kepler_Tcent_T0P(
                sample_plan[:, convert_out['P']], sample_plan[:, convert_out['f']],
                sample_plan[:, convert_out['e']], sample_plan[:, convert_out['o']])

            if dynamical_flag:
                convert_out['K'] = n_orbital + 5
                sample_plan[:, convert_out['K']] = kp.kepler_K1(mc.star_mass[0],
                                                                sample_plan[:, convert_out['M']]/mc.M_SEratio,
                                                                sample_plan[:, convert_out['P']],
                                                                sample_plan[:, convert_out['i']],
                                                                sample_plan[:, convert_out['e']])

            if 'curvature' in mc.model_list:
                for ii in xrange(0, nsample):
                    convert_tmp = mc.models[model].convert(flatchain[ii, :])
                    for var in mc.models[model].list_pams:
                        sample_plan[ii, convert_out[var]] = convert_tmp[var]

                #for var in mc.models[model].list_pams:
                #    sample_med[convert_out[var], 0] = np.median(sample_plan[:, convert_out[var]])

            for n, (P, K, e, M) in enumerate(zip(sample_plan[:, convert_out['P']],
                                                 sample_plan[:, convert_out['K']],
                                                 sample_plan[:, convert_out['e']],
                                                 M_star1_rand)):
                # Planet mass
                sample_plan[n, convert_out['M_kep']] = mc.M_SJratio * scipy.optimize.fsolve(get_mass, x0, args=(M, P, K, e))
                # semi-major axis
                sample_plan[n, convert_out['a_smj']] = np.power(
                    (Mu_sun * np.power(P * seconds_in_day / (2 * np.pi), 2) / (AU_km ** 3.0)) * M, 1.00 / 3.00)

            if not dynamical_flag:
                if mc.models['planets'].inclination[planet_name][1] > 0.01:
                    sample_plan[:, convert_out['M_kep']] = sample_plan[:, convert_out['M_kep']] / np.sin(np.pi / 180. *
                            np.random.normal(mc.models['planets'].inclination[planet_name][0], mc.models['planets'].inclination[planet_name][1], nsample))

            """ Rescale the center point and range for omega and phase"""
            for key in ['o','f']:
                var_index = convert_out[key]
                sel_turnA = (sample_plan[:, var_index] > np.pi + median_tmp[var_index])
                sel_turnB = (sample_plan[:, var_index] < median_tmp[var_index] - np.pi)
                sample_plan[sel_turnA, var_index] -= 2*np.pi
                sample_plan[sel_turnB, var_index] += 2*np.pi

            sample_med = np.asarray(map(lambda v: (v[1], v[2] - v[1], v[1] - v[0]),
                                        zip(*np.percentile(sample_plan[:, :], [15.865, 50, 84.135], axis=0))))

            #sample_med[o_index, 0] = sample_tmp[o_index,0]
            #sample_plan_tmp = sample_plan[o_index, :]

            #sample_med[o_index, 1:] = np.abs(np.percentile(sample_plan[o_index, :], [15.865, 84.135]))-sample_med[o_index, 0]

            e_med = np.percentile(sample_plan[:, convert_out['e']], 65.865, axis=0)

            if sampler in sample_keyword['polychord']:
                ''' Let's save all the samples for this specific planet, to be used later in the special PolyChrod plots'''
                sample_total[planet_name] = {'sample_plan': sample_plan, 'convert_out': convert_out}

            print 'Period = ', sample_med[convert_out['P'], 0], ' +\sigma ', sample_med[convert_out['P'], 1], ' -\sigma ', sample_med[convert_out['P'], 2]
            print 'K      = ', sample_med[convert_out['K'], 0], ' +\sigma ', sample_med[convert_out['K'], 1], ' -\sigma ', sample_med[convert_out['K'], 2]
            print 'phase  = ', sample_med[convert_out['f'], 0], ' +\sigma ', sample_med[convert_out['f'], 1], ' -\sigma ', sample_med[convert_out['f'], 2]
            print 'e      = ', sample_med[convert_out['e'], 0], ' +\sigma ', sample_med[convert_out['e'], 1], ' -\sigma ', sample_med[convert_out['e'], 2], ', < ', e_med
            print 'o      = ', sample_med[convert_out['o'], 0], ' +\sigma ', sample_med[convert_out['o'], 1], ' -\sigma ', sample_med[convert_out['o'], 2]

            if dynamical_flag:
                print 'lN     = ', sample_med[convert_out['lN'], 0], ' +\sigma ', sample_med[convert_out['lN'], 1], ' -\sigma ', sample_med[convert_out['lN'], 2]
                print 'i      = ', sample_med[convert_out['i'], 0], ' +\sigma ', sample_med[convert_out['i'], 1], ' -\sigma ', sample_med[convert_out['i'], 2]
                print 'Mass_J = ', sample_med[convert_out['M'], 0]/mc.M_JEratio, ' +\sigma ', sample_med[convert_out['M'], 1]/mc.M_JEratio, ' -\sigma ', sample_med[convert_out['M'], 2]/mc.M_JEratio
                print 'Mass_E = ', sample_med[convert_out['M'], 0], ' +\sigma ', sample_med[convert_out['M'], 1], ' -\sigma ', sample_med[convert_out['M'], 2]
            else:
                print 'Mass_J = ', sample_med[convert_out['M_kep'], 0], ' +\sigma ', sample_med[convert_out['M_kep'], 1], ' -\sigma ', sample_med[convert_out['M_kep'], 2]
                print 'Mass_E = ', sample_med[convert_out['M_kep'], 0]*mc.M_JEratio, ' +\sigma ',sample_med[convert_out['M_kep'], 1]*mc.M_JEratio, ' -\sigma ',sample_med[convert_out['M_kep'], 2]*mc.M_JEratio

            print 'Tperi  = ', sample_med[convert_out['Tperi'], 0], ' +\sigma ', sample_med[convert_out['Tperi'], 1], ' -\sigma ', sample_med[convert_out['Tperi'], 2]
            print 'Tcent  = ', sample_med[convert_out['Tcent'], 0], ' +\sigma ', sample_med[convert_out['Tcent'], 1], ' -\sigma ', sample_med[convert_out['Tcent'], 2]
            print 'a      = ', sample_med[convert_out['a_smj'], 0], ' +\sigma ', sample_med[convert_out['a_smj'], 1], ' -\sigma ', sample_med[convert_out['a_smj'], 2]
            print

            sel_list = []
            sel_label = []

            if 'P' in mc.variable_list[planet_name]:
                sel_list.append(convert_out['P'])
                sel_label.append('P [d]')
            if 'K' in mc.variable_list[planet_name]:
                sel_list.append(convert_out['K'])
                sel_label.append('K [$m s_{-1}$]')

            if 'f' in mc.variable_list[planet_name]:
                sel_list.append(convert_out['f'])
                sel_label.append('$\phi$ [rad]')

            if 'e' in mc.variable_list[planet_name]:
                sel_list.append(convert_out['e'])
                sel_label.append('e')
            if 'o' in mc.variable_list[planet_name]:
                sel_list.append(convert_out['o'])
                sel_label.append('$\omega$ [rad]')
            if 'ecoso' in mc.variable_list[planet_name]:
                sel_list.append(convert_out['e'])
                sel_label.append('e')
                sel_list.append(convert_out['o'])
                sel_label.append('$\omega$ [rad]')

            if 'M' in mc.variable_list[planet_name]:
                sel_list.append(convert_out['M'])
                sel_label.append('M [$M_\oplus $')
            else:
                sel_list.append(convert_out['M_kep'])
                #sample_plan[:,convert_out['M_kep']] *= mc.M_JEratio
                sel_label.append('M [$M_j$]')

            if 'curvature' in mc.model_list:
                for n_var, var in enumerate(mc.models[model].list_pams):
                    sel_list.append(convert_out[var])
                    sel_label.append('C$_'+repr(n_var)+'$')

            fig = corner.corner(sample_plan[:, sel_list], labels=sel_label, truths=sample_med[sel_list, 0])
            fig.savefig(dir_output + planet_name + "_corners.pdf", bbox_inches='tight', dpi=300)
            plt.close(fig)

            if args.p != 'False' or args.v != 'False':
                # Write down the residuals
                color_list = ['b', 'g', 'r', 'c', 'm', 'y', 'k']

                f1, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=True, figsize=(12, 12))
                f2, (ax3, ax4) = plt.subplots(2, sharex=True, sharey=True, figsize=(12, 12))

                ax1.set_xlabel('BJD [d]')
                ax3.set_xlabel('Orbital phase')

                ax1.set_ylabel('RV [$m s_{-1}$]')
                ax2.set_ylabel('RV residuals [$m s_{-1}$]')
                ax3.set_ylabel('RV [$m s_{-1}$]')
                ax4.set_ylabel('RV residuals [$m s_{-1}$]')

                ax1.plot(x_range, model_plan['BJD'][planet_name], c='g')
                ax3.plot(x_phase, model_plan['pha'][planet_name], c='g')

                color_count = 0
                for dataset_name, dataset in mc.dataset_dict.items():
                    if dataset.kind == 'RV':
                        col_sel = color_list[color_count % 7]
                        color_count += 1
                        p_pha = (dataset.x0 / sample_med[convert_out['P'], 0]) % 1
                        y_det = dataset.y - model_dsys[dataset_name] - model_curv[dataset_name]
                        y_res = dataset.y - model_dsys[dataset_name] - \
                                model_orbs[dataset_name] - model_curv[dataset_name] -\
                                model_actv[dataset_name]

                        y_1pl = y_res + model_plan[dataset_name][planet_name]

                        ax1.errorbar(dataset.x, y_1pl, yerr=dataset.e, fmt=col_sel + '.', zorder=2)
                        ax2.errorbar(dataset.x, y_res, yerr=dataset.e, fmt=col_sel + '.', zorder=2)

                        ax3.errorbar(p_pha, y_1pl, yerr=dataset.e, fmt=col_sel + '.', zorder=2)
                        ax4.errorbar(p_pha, y_res, yerr=dataset.e, fmt=col_sel + '.', zorder=2)

                        fileout = open(plot_dir + planet_name + '_' + dataset_name + '_kep.dat', 'w')
                        fileout.write('descriptor BJD pha RV,+- RVdet,+- RVpla,+- RVres,+- RVmod,+- \n')
                        for ii in xrange(0, dataset.n):
                            fileout.write('{0:14f} {1:14f} {2:14f} {3:14f} {4:14f} {5:14f} {6:14f} {7:14f} '
                                          '{8:14f} {9:14f} {10:14f} {11:14f}'
                                          '\n'.format(dataset.x[ii], p_pha[ii],
                                                      dataset.y[ii], dataset.e[ii],
                                                      y_det[ii], dataset.e[ii], y_1pl[ii], dataset.e[ii],
                                                      y_res[ii], dataset.e[ii],
                                                      model_orbs[dataset_name][ii], dataset.e[ii]))
                        fileout.close()

                f1.subplots_adjust(hspace=0)
                f1.savefig(plot_dir + planet_name + '_kep.pdf', bbox_inches='tight', dpi=300)
                f2.subplots_adjust(hspace=0)
                f2.savefig(plot_dir + planet_name + '_pha.pdf', bbox_inches='tight', dpi=300)
                plt.close(f1)
                plt.close(f2)

                fileout = open(plot_dir + planet_name + '_kep.dat', 'w')
                fileout.write('descriptor x_range m_kepler \n')
                for ii in xrange(0, np.size(x_range)):
                    fileout.write('{0:14f} {1:14f} \n'.format(x_range[ii], model_plan['BJD'][planet_name][ii]))
                fileout.close()

                fileout = open(plot_dir + planet_name + '_pha.dat', 'w')
                fileout.write('descriptor x_phase m_phase \n')
                for ii in xrange(0, np.size(x_phase)):
                    fileout.write('{0:14f} {1:14f} \n'.format(x_phase[ii], model_plan['pha'][planet_name][ii]))
                fileout.close()

                if args.v != 'False':

                    veusz_dir = dir_output + '/Veuz_plot/'
                    print veusz_dir
                    if not os.path.exists(veusz_dir):
                        os.makedirs(veusz_dir)

                    list_labels = ['P', 'K', 'e', 'M_kep']

                    n_int = 4

                    output_plan = np.zeros([nsample, n_int], dtype=np.double)
                    output_plan[:, 0] = sample_plan[:, convert_out['P']]
                    output_plan[:, 1] = sample_plan[:, convert_out['K']]
                    output_plan[:, 2] = sample_plan[:, convert_out['e']]
                    output_plan[:, 3] = sample_plan[:, convert_out['M_kep']]
                    plot_truths = np.percentile(output_plan[:, :], [15.865, 50, 84.135], axis=0)

                    #veusz_dir = dir_output + '/Veuz_plot/'
                    #if not os.path.exists(veusz_dir):
                    #    os.makedirs(veusz_dir)

                    n_bins = 60 + 1

                    h5f = h5py.File(veusz_dir + planet_name + '_hist1d.hdf5', "w")
                    data_grp = h5f.create_group("hist1d")

                    data_lim = np.zeros([n_int, 2], dtype=np.double)
                    data_edg = np.zeros([n_int, n_bins], dtype=np.double)
                    for ii in xrange(0, n_int):
                        data_lim[ii, :] = [np.amin(output_plan[:, ii]), np.amax(output_plan[:, ii])]
                        if data_lim[ii, 0] == data_lim[ii, 1]:
                            data_lim[ii, :] = [0.0, 0.1]
                        data_edg[ii, :] = np.linspace(data_lim[ii, 0], data_lim[ii, 1], n_bins)

                    for ii in xrange(0, n_int):
                        for jj in xrange(ii, n_int):
                            x_data = output_plan[:, ii]
                            y_data = output_plan[:, jj]
                            x_edges = data_edg[ii, :]
                            y_edges = data_edg[jj, :]

                            if ii != jj:
                                hist2d = np.histogram2d(x_data, y_data, bins=[x_edges, y_edges], normed=True)
                                hist1d_y = np.histogram(y_data, bins=y_edges, normed=True)

                                Hflat = hist2d[0].flatten()
                                inds = np.argsort(Hflat)[::-1]
                                Hflat = Hflat[inds]
                                sm = np.cumsum(Hflat)
                                sm /= sm[-1]

                                x_edges_1d = (x_edges[1:] + x_edges[:-1])/2
                                y_edges_1d = (y_edges[1:] + y_edges[:-1])/2
                                h2d_out = np.zeros([n_bins, n_bins])
                                h2d_out[0, 1:] = x_edges_1d
                                h2d_out[1:, 0] = y_edges_1d
                                h2d_out[1:, 1:] = hist2d[0].T *1. / np.amax(hist2d[0])

                                h2d_list =  h2d_out.tolist()
                                h2d_list[0][0] = ''
                                csvfile = veusz_dir + planet_name + '_hist2d_' + list_labels[ii] + '_' + list_labels[jj] + '.csv'
                                with open(csvfile, "w") as output:
                                    writer = csv.writer(output, lineterminator='\n')
                                    writer.writerows(h2d_list)
                            else:
                                hist1d = np.histogram(x_data, bins=x_edges)
                                hist1d_norm = hist1d[0]*1. / nsample
                                x_edges_1d = (x_edges[1:]+ x_edges[:-1])/2
                                data_grp.create_dataset(list_labels[ii]+'_x', data=x_edges_1d, compression="gzip")
                                data_grp.create_dataset(list_labels[ii]+'_y', data=hist1d_norm, compression="gzip")

                    # plot lower-upper limits for the mass

                    pams_limits = np.zeros([n_total, 2])
                    pams_limits[:, 0] = sample_med[:, 0] - sample_med[:, 1]
                    pams_limits[:, 1] = sample_med[:, 0] + sample_med[:, 2]

                    if pams_limits[convert_out['e'], 0] < 0.02 : pams_limits[convert_out['e'], 0] = 0.00
                    y_flg = np.ones(nsample, dtype=bool)

                    for var in ['P', 'K', 'e', 'f', 'o']:
                        y_flg &= (sample_plan[:, convert_out[var]] >= pams_limits[convert_out[var], 0]) & \
                                (sample_plan[:, convert_out[var]] <= pams_limits[convert_out[var], 1])
                        print ' LIMITS ', var, pams_limits[convert_out[var], :]

                    y_flg &= (flatlnprob[:] > lnprob_med[0] - lnprob_med[1]) # & (flatlnprob[:] < lnprob_med[0] + lnprob_med[2])

                    random_n = 10000
                    if np.sum(y_flg) <= random_n:
                        random_n = np.sum(y_flg)
                        ii_kep = np.where(y_flg)[0][:]
                    else:
                        ii_kep = np.random.permutation(np.where(y_flg)[0][:])[:random_n]
                        #print np.where(y_flg)
                        #ii_selected = np.where(y_flg)
                        #print ii_selected
                        #ii_kep = ii_selected[np.random.randint(low=0, high=np.sum(y_flg), size=random_n).tolist()]

                    #ii_kep = np.random.randint(low=0, high=nsample, size=random_n)

                    x_kep = np.arange(5900, 7900, 2)
                    y_kep = np.zeros([np.size(x_kep), 3])

                    x_pha = np.arange(-1.00, 2.00, 0.005,dtype=np.double)
                    y_pha = np.zeros([np.size(x_pha), 3])

                    #y_flg = np.ones(random_n, dtype=bool)

                    y_kep_tmp = np.zeros(np.size(x_kep))
                    y_pha_tmp = np.zeros(np.size(x_pha))

                    #if pams_limits[convert_out['e'], 0] < 0.02 : pams_limits[convert_out['e'], 0] = 0.00
                    #for var in convert_out:
                    #    y_flg = y_flg & (sample_plan[ii_kep, convert_out[var]] >= pams_limits[convert_out[var], 0]) & \
                    #            (sample_plan[ii_kep, convert_out[var]] <= pams_limits[convert_out[var], 1])
                    #    print ' LIMITS ', pams_limits[convert_out[var], :]
                    #y_flg = y_flg & (flatlnprob[ii_kep] > lnprob_med[0] - lnprob_med[1]) & (flatlnprob[ii_kep] < lnprob_med[0] + lnprob_med[2])

                    y_kep[:, 2] = kp.kepler_RV_T0P(x_kep - mc.Tref,
                                                   sample_med[convert_out['f'], 0],
                                                   sample_med[convert_out['P'], 0],
                                                   sample_med[convert_out['K'], 0],
                                                   sample_med[convert_out['e'], 0],
                                                   sample_med[convert_out['o'], 0])
                    y_kep[:, 0] = y_kep[:, 2]
                    y_kep[:, 1] = y_kep[:, 2]
                    # value initialization
                    #y_pha[:, 2] = kp.kepler_RV_T0P(x_pha*sample_med[0, 0], sample_med[2, 0], sample_med[0, 0], sample_med[1, 0],
                    #                               sample_med[3, 0], sample_med[4, 0])

                    y_pha[:, 2] = kp.kepler_RV_T0P(x_pha*sample_med[convert_out['P'], 0],
                                                   sample_med[convert_out['f'], 0],
                                                   sample_med[convert_out['P'], 0],
                                                   sample_med[convert_out['K'], 0],
                                                   sample_med[convert_out['e'], 0],
                                                   sample_med[convert_out['o'], 0])
                    y_pha[:, 0] = y_pha[:, 2]
                    y_pha[:, 1] = y_pha[:, 2]

                    print 'Accepted randomizations: ', np.sum(y_flg)
                    print 'Kept randomizations: ', random_n

                    for ir in ii_kep:
                        #if y_flg[ii]:
                            #y_kep_tmp = kp.kepler_RV_T0P(x_kep - mc.Tref, sample_plan[ir, 2], sample_plan[ir, 0], sample_plan[ir, 1], sample_plan[ir, 3], sample_plan[ir, 4])
                            #y_pha_tmp = kp.kepler_RV_T0P(x_pha*sample_plan[ir, 0], sample_plan[ir, 2], sample_plan[ir, 0], sample_plan[ir, 1], sample_plan[ir, 3], sample_plan[ir, 4])

                            y_kep_tmp = kp.kepler_RV_T0P(x_kep - mc.Tref,
                                                         sample_plan[ir, convert_out['f']],
                                                         sample_plan[ir, convert_out['P']],
                                                         sample_plan[ir, convert_out['K']],
                                                         sample_plan[ir, convert_out['e']],
                                                         sample_plan[ir, convert_out['o']])
                            y_pha_tmp = kp.kepler_RV_T0P(x_pha*sample_plan[ir, convert_out['P']],
                                                         sample_plan[ir, convert_out['f']],
                                                         sample_plan[ir, convert_out['P']],
                                                         sample_plan[ir, convert_out['K']],
                                                         sample_plan[ir, convert_out['e']],
                                                         sample_plan[ir, convert_out['o']])

                            y_kep[:, 0] = np.minimum(y_kep[:, 0], y_kep_tmp)
                            y_kep[:, 1] = np.maximum(y_kep[:, 1], y_kep_tmp)

                            y_pha[:, 0] = np.minimum(y_pha[:, 0], y_pha_tmp)
                            y_pha[:, 1] = np.maximum(y_pha[:, 1], y_pha_tmp)

                    fig0 = plt.figure(0, figsize=(12, 12))
                    plt.figure(0)
                    plt.plot(x_kep, y_kep[:, 2], c='g')
                    plt.plot(x_kep, y_kep[:, 0], c='r')
                    plt.plot(x_kep, y_kep[:, 1], c='b')
                    plt.savefig(veusz_dir + planet_name + '_kep.png', bbox_inches='tight', dpi=300)
                    plt.close(fig0)

                    fig1 = plt.figure(1, figsize=(12, 12))
                    plt.figure(1)
                    plt.plot(x_pha, y_pha[:, 2], c='g')
                    plt.plot(x_pha, y_pha[:, 0], c='r')
                    plt.plot(x_pha, y_pha[:, 1], c='b')
                    plt.savefig(veusz_dir + planet_name + '_pha.png', bbox_inches='tight', dpi=300)
                    plt.close(fig1)

                    # h5f = h5py.File('output/'+planet_name+"planet"+`pp`+'_kep.hdf5', "w")
                    # data_grp = h5f.create_group("data")
                    # data_grp.create_dataset("x_kep", data=x_kep, compression="gzip")
                    # data_grp.create_dataset("y_kep", data=y_kep, compression="gzip")
                    # h5f.close()
                    # h5f = h5py.File('output/'+planet_name+"planet"+`pp`+'_pha.hdf5', "w")
                    # data_grp = h5f.create_group("data")
                    # data_grp.create_dataset("x_pha", data=x_kep, compression="gzip")
                    # data_grp.create_dataset("y_pha", data=y_kep, compression="gzip")
                    # h5f.close()

                    fileout = open(veusz_dir + planet_name + '_kep.dat','w')
                    fileout.write('descriptor x_kep m_kep y_kep,+- \n')
                    for ii in xrange(0,np.size(x_kep)):
                        fileout.write('{0:14f} {1:14f} {2:14f} {3:14f}  \n'.format(\
                            x_kep[ii], y_kep[ii, 2], (y_kep[ii, 1]+y_kep[ii, 0])/2, (y_kep[ii, 1]-y_kep[ii, 0])/2))
                    fileout.close()

                    fileout = open(veusz_dir + planet_name + '_pha.dat','w')
                    fileout.write('descriptor x_pha m_pha y_pha,+- \n')
                    for ii in xrange(0,np.size(x_pha)):
                        fileout.write('{0:14f} {1:14f} {2:14f} {3:14f} \n'.format(\
                            x_pha[ii], y_pha[ii, 2], (y_pha[ii, 1]+y_pha[ii, 0])/2, (y_pha[ii, 1]-y_pha[ii, 0])/2))
                    fileout.close()

            print 'Planet ', planet_name, ' completed'
            print
            print '-----------------------------'
            print

        if sampler in sample_keyword['polychord']:
            ''' Now we do the Polychord plots were all the points are collected together'''

            color_list = ['b', 'c', 'm', 'g', 'r', 'k', 'y']

            i_color = 0

            fig0 = plt.figure(0, figsize=(12, 12))
            fig1 = plt.figure(1, figsize=(12, 12))
            fig2 = plt.figure(2, figsize=(12, 12))
            fig3 = plt.figure(3, figsize=(12, 12))
            fig4 = plt.figure(4, figsize=(12, 12))

            for planet_name in mc.models['planets'].planet_name:

                dynamical_flag = (planet_name in mc.models['planets'].dynamical)

                convert_out = sample_total[planet_name]['convert_out']
                if dynamical_flag:
                    sample_M = sample_total[planet_name]['sample_plan'][:, convert_out['M']]
                else:
                    sample_M = sample_total[planet_name]['sample_plan'][:, convert_out['M_kep']]*mc.M_JEratio

                sample_P = sample_total[planet_name]['sample_plan'][:, convert_out['P']]
                sample_e = sample_total[planet_name]['sample_plan'][:, convert_out['e']]

                plt.figure(0)
                plt.scatter(sample_P, sample_M, color=color_list[i_color], alpha=0.2, s=5, linewidths=None)
                plt.draw()
                plt.figure(1)
                plt.scatter(sample_M, sample_e, color=color_list[i_color], alpha=0.2, s=5, linewidths=None)
                plt.draw()
                plt.figure(2)
                plt.scatter(sample_P, sample_e, color=color_list[i_color], alpha=0.2, s=5, linewidths=None)
                plt.draw()
                plt.figure(3)
                plt.scatter(sample_P, flatlnprob, color=color_list[i_color], alpha=0.2, s=5, linewidths=None)
                plt.draw()
                plt.figure(4)
                for planet_alt in mc.models['planets'].planet_name:
                    if planet_alt == planet_name: continue
                    sample_P_alt = sample_total[planet_alt]['sample_plan'][:, convert_out['P']]
                    plt.scatter(sample_P, sample_P_alt, color=color_list[i_color], alpha=0.2, s=5, linewidths=None)
                plt.draw()

                i_color += 1
                if i_color == np.size(color_list):
                    i_color = 0

            plt.figure(0)
            plt.xlabel('P [d]')
            plt.ylabel('M [$M_\oplus $]')
            plt.xscale("log", nonposx='clip')
            plt.draw()
            plt.savefig(dir_output + 'PolyScatter_P_M.pdf', bbox_inches='tight', dpi=300)
            plt.close(fig0)
            plt.figure(1)
            plt.xlabel('M [$M_\oplus $]')
            plt.ylabel('e')
            plt.draw()
            plt.savefig(dir_output + 'PolyScatter_M_e.pdf', bbox_inches='tight', dpi=300)
            plt.close(fig1)
            plt.figure(2)
            plt.xlabel('P [d]')
            plt.ylabel('e')
            plt.xscale("log", nonposx='clip')
            plt.draw()
            plt.savefig(dir_output + 'PolyScatter_P_e.pdf', bbox_inches='tight', dpi=300)
            plt.close(fig2)
            plt.figure(3)
            plt.xlabel('P [d]')
            plt.ylabel('lnP')
            plt.xscale("log", nonposx='clip')
            plt.draw()
            plt.savefig(dir_output + 'PolyScatter_P_lnP.pdf', bbox_inches='tight', dpi=300)
            plt.close(fig3)
            plt.figure(4)
            plt.xlabel('P [d]')
            plt.ylabel('e')
            plt.xscale("log", nonposx='clip')
            plt.yscale("log", nonposx='clip')
            plt.draw()
            plt.savefig(dir_output + 'PolyScatter_P_P.pdf', bbox_inches='tight', dpi=300)
            plt.close(fig4)

        if args.forecast is not None:
            y_flg = (flatlnprob[:] > lnprob_med[0] - lnprob_med[1]) # & (flatlnprob[:] < lnprob_med[0] + lnprob_med[2])

            random_n = 1000
            if np.sum(y_flg) <= random_n:
                random_n = np.sum(y_flg)
                ii_kep = np.where(y_flg)
            else:
                ii_kep = np.random.permutation(np.where(y_flg)[0][:])[:random_n]

            print 'Accepted randomizations: ', np.sum(y_flg)
            print 'Kept randomizations: ', random_n
            print np.size(flatchain[0, :]), np.size(flatchain[:, 0])

            x_kep = np.arange(boundaries[0], args.forecast, 1.)
            x_pha = np.arange(0.0, 1.0, 0.01)
            y_kep = np.zeros([np.size(x_kep), 3])

            for ir in ii_kep:
                random_dsys, random_plan, random_orbs, random_actv, random_curv = mc.rv_make_model(flatchain[ir, :], x_kep, x_pha)
                model_total = random_dsys['BJD'] + random_actv['BJD'] + random_orbs['BJD'] + random_curv['BJD']
                plt.plot(x_kep, model_total, c='k', alpha=0.05)
            plt.savefig(dir_output+'_forecast.pdf', bbox_inches='tight', dpi=300)
            plt.close()

    if mc.models[model].model_class == 'gaussian_process':

        print 'Gaussian process summary'
        print

        n_vars = 0
        sample_plan_transpose = []
        sel_label = []

        for name in mc.models[model].list_pams_common:
            if name in mc.variable_list[mc.models[model].common_ref]:
                n_vars += 1
                # mc.models[model].var_list[name] select the corresponding value in the emcee input array (theta)
                var = flatchain[:, mc.models[model].var_list[mc.models[model].common_ref][name]]
                # mc.models[model].variables[name](theta, mc.models[model].fixed, :) convert the value into physical unit
                # e.g. from logarithmic to linear space)
                var_phys = mc.models[model].variables[mc.models[model].common_ref][name](var, var, xrange(0, nsample))
                sample_plan_transpose.append(var_phys)
                sel_label.append(name)

        for dataset_name, dataset in mc.dataset_dict.items():
            for model in dataset.models:
                if mc.models[model].model_class == 'gaussian_process':
                    for name in mc.models[model].list_pams_dataset:
                        n_vars += 1
                        var = flatchain[:, mc.models[model].var_list[dataset_name][name]]
                        var_phys = mc.models[model].variables[dataset_name][name](var, var, xrange(0, nsample))
                        sample_plan_transpose.append(var_phys)
                        sel_label.append('$'+ dataset.name + '$ ' + name)

        sample_plan = np.asarray(sample_plan_transpose).T
        sample_med = np.asarray(map(lambda v: (v[1], v[2] - v[1], v[1] - v[0]),
                                        zip(*np.percentile(sample_plan[:, :], [15.865, 50, 84.135], axis=0))))

        for lab, sam in zip(sel_label, sample_med):
            print lab,' = ', sam[0], ' +\sigma ', sam[1], ' -\sigma ', sam[2]
        print

        fig = corner.corner(sample_plan[:, :], labels=sel_label, truths=sample_med[:, 0])
        fig.savefig(dir_output + "GPs_corners.pdf", bbox_inches='tight', dpi=300)
        plt.close(fig)

        print 'Gaussian process summary completed'
        print
        print '-----------------------------'
        print

    if mc.models[model].model_class == 'curvature':

        print 'Curvature summary'
        print

        n_vars = 0
        sample_plan_transpose = []
        sel_label = []

        for name in mc.models[model].list_pams:
            n_vars += 1
            var = flatchain[:, mc.models[model].var_list[name]]
            var_phys = mc.models[model].variables[name](var, var, xrange(0, nsample))
            sample_plan_transpose.append(var_phys)
            sel_label.append(name)

        sample_plan = np.asarray(sample_plan_transpose).T
        sample_med = np.asarray(map(lambda v: (v[1], v[2] - v[1], v[1] - v[0]),
                                        zip(*np.percentile(sample_plan[:, :], [15.865, 50, 84.135], axis=0))))

        for lab, sam in zip(sel_label, sample_med):
            print lab,' = ', sam[0], ' +\sigma ', sam[1], ' -\sigma ', sam[2]
        print

        fileout = open(plot_dir + 'curvature.dat', 'w')
        fileout.write('descriptor x_range m_curv \n')
        for ii in xrange(0, np.size(x_range)):
            fileout.write('{0:14f} {1:14f} \n'.format(x_range[ii], model_curv['BJD'][ii]))
        fileout.close()

        print 'Curvature summary completed'
        print
        print '-----------------------------'
        print

    if mc.models[model].model_class == 'curvature':

        print 'Curvature summary'
        print

        n_vars = 0
        sample_plan_transpose = []
        sel_label = []

        for name in mc.models[model].list_pams:
            n_vars += 1
            var = flatchain[:, mc.models[model].var_list[name]]
            var_phys = mc.models[model].variables[name](var, var, xrange(0, nsample))
            sample_plan_transpose.append(var_phys)
            sel_label.append(name)

        sample_plan = np.asarray(sample_plan_transpose).T
        sample_med = np.asarray(map(lambda v: (v[1], v[2] - v[1], v[1] - v[0]),
                                    zip(*np.percentile(sample_plan[:, :], [15.865, 50, 84.135], axis=0))))

        for lab, sam in zip(sel_label, sample_med):
            print lab, ' = ', sam[0], ' +\sigma ', sam[1], ' -\sigma ', sam[2]
        print

        fileout = open(plot_dir + 'curvature.dat', 'w')
        fileout.write('descriptor x_range m_curv \n')
        for ii in xrange(0, np.size(x_range)):
            fileout.write('{0:14f} {1:14f} \n'.format(x_range[ii], model_curv['BJD'][ii]))
        fileout.close()

        print 'Curvature summary completed'
        print
        print '-----------------------------'
        print

        print 'Correlation summary'
        print

        n_vars = 0
        sample_plan_transpose = []
        sel_label = []
        plt.rc('text', usetex=False)

        for name_ref in mc.models[model].list_pams:
            for name_asc in mc.models[model].list_pams[name_ref]:
                for name_var in mc.models[model].list_pams[name_ref][name_asc]:
                    n_vars += 1
                    var = flatchain[:, mc.models[model].var_list[name_ref][name_asc][name_var]]
                    var_phys = mc.models[model].variables[name_ref][name_asc][name_var](var, var, xrange(0, nsample))
                    sample_plan_transpose.append(var_phys)
                    sel_label.append(name_ref+'_'+name_asc+'_'+name_var)

        sample_plan = np.asarray(sample_plan_transpose).T
        sample_med = np.asarray(map(lambda v: (v[1], v[2] - v[1], v[1] - v[0]),
                                        zip(*np.percentile(sample_plan[:, :], [15.865, 50, 84.135], axis=0))))

        for name_ref in mc.models[model].list_pams:
            for name_asc in mc.models[model].list_pams[name_ref]:
                print name_ref, name_asc, 'zero-point', mc.models[model].x_zero[name_ref][name_asc]

        for lab, sam in zip(sel_label, sample_med):
            print lab,' = ', sam[0], ' +\sigma ', sam[1], ' -\sigma ', sam[2]
        print



        fileout = open(plot_dir + 'correlation.dat', 'w')
        fileout.write('descriptor x_range m_corr \n')
        #for ii in xrange(0, np.size(x_range)):
        #    fileout.write('{0:14f} {1:14f} \n'.format(x_range[ii], model_act['BJD'][ii]))
        #fileout.close()

        fig = corner.corner(sample_plan[:, :], labels=sel_label, truths=sample_med[:, 0])
        fig.savefig(dir_output + "Correlation_corners.pdf", bbox_inches='tight', dpi=300)
        plt.close(fig)

        plt.rc('text', usetex=True)

        print 'Curvature summary completed'
        print
        print '-----------------------------'
        print

print
