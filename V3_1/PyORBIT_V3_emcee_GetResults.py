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
import sys
mpl.use('Agg')
from matplotlib import pyplot as plt
import corner
sys.path.append('/Users/malavolta/Astro/CODE/trades/pytrades')
#from pytrades_lib import pytrades
import constants

def GelmanRubin_old(chains_input):
    n_iters, n_chain = np.shape(chains_input)
    W = np.asarray(0., dtype=np.double)
    z_pc = np.sum(chains_input, axis=0) / n_iters  # eq 20
    for nn in xrange(0, n_chain):
        W += (np.sum(np.power(chains_input[:, nn] - z_pc[nn], 2))) / ((n_iters - 1) * n_chain)  # eq 21
    z_pp = np.sum(chains_input) / (n_chain * n_iters)
    B = np.sum(np.power(z_pc - z_pp, 2)) * (n_chain / (n_iters - 1))
    var = W * (n_chain - 1) / n_chain + B / n_chain
    return np.sqrt(var / W)

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


parser = argparse.ArgumentParser(prog='PyORBIT_V3_GetResults.py', description='Extract results from output MCMC')
# parser.add_argument('-l', type=str, nargs='+', help='line identificator')
parser.add_argument('config_file', type=str, nargs=1, help='config file')
parser.add_argument('-p', type=str, nargs='?', default='False', help='Create plot files')
parser.add_argument('-mp', type=str, nargs='?', default='False', help='Create MEGA plot')
parser.add_argument('-v', type=str, nargs='?', default='False', help='Create Veusz ancillary files')
parser.add_argument('-t', type=str, nargs='?', default='False', help='Create GR traces')
parser.add_argument('-nburn', type=int, nargs='?', default=0, help='emcee burn-ins')
parser.add_argument('-cc', type=str, nargs='?', default='False', help='Use ChainConsumer')

args = parser.parse_args()

file_conf = args.config_file[0]

# file_conf = raw_input()

mc = ModelContainer()
yaml_parser(file_conf, mc)

mc.create_bounds()
if bool(mc.pcv.dynamical):
    mc.pcv.prepare_dynamical(mc)

M_star1 = mc.star_mass[0]
M_star1_err = mc.star_mass[1]

dir_output = './' + mc.planet_name + '/'

mc.variable_list = pickle.load(open(dir_output + 'vlist.pick', 'rb'))
mc.scv.use_offset = pickle.load(open(dir_output + 'scv_offset.pick', 'rb'))

print mc.variable_list

h5f = h5py.File(dir_output + mc.planet_name + '.hdf5', "r")

h5f_data = h5f['/data']
h5f_emcee = h5f['/emcee']

for item in h5f_emcee.attrs.keys():
    print item + ":", h5f_emcee.attrs[item]

    if item == 'nwalkers': mc.nwalkers = h5f_emcee.attrs[item]
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
#print 'Autocorrelation time'
#print acor[:]
print
print '*************************************************************'
print

if args.nburn > 0:
    mc.nburn = args.nburn

if mc.nsave >0 :
    mc.nsteps = h5f_emcee.attrs['nsample']
    if mc.nburn > mc.nsteps:
        mc.nburn = mc.nsteps / 4

mc.nsample = mc.nsteps / mc.thin
mc.nburnin = mc.nburn / mc.thin

lnprb_T = lnprobability[:][:].T

chain_T = np.ndarray([mc.nsample, mc.nwalkers, mc.ndim], dtype=np.double)
for ii in xrange(0, mc.ndim):
    chain_T[:, :, ii] = chain[:, :, ii].T

chain_burnt = chain_T[mc.nburnin:, :, :]
s = chain_burnt.shape
lnprob_burnt = lnprb_T[mc.nburnin:, :,]
flatchain = chain_burnt.reshape(s[1] * s[0], s[2])
flatlnprob = lnprob_burnt.reshape(s[1] * s[0])
n_kept = s[1] * s[0]
#sel_flatchain = flatchain[:, 0] < 1.

chain_med = np.asarray(map(lambda v: (v[1], v[2] - v[1], v[1] - v[0]),
                           zip(*np.percentile(flatchain[:, :], [15.865, 50, 84.135], axis=0))))

mc.results_resumen(chain_med[:, 0])

print
print '*************************************************************'
print

lnprob_median = np.median(flatlnprob)
fig = plt.plot(lnprb_T[:, :], '-', alpha=0.5)
plt.axhline(lnprob_median)
plt.axvline(mc.nburnin, c='r')
plt.savefig(dir_output + 'LNprob_chain.png', bbox_inches='tight')
plt.close()
print 'LNprob median = ', lnprob_median

print
print 'Tref: ', mc.Tref
print
print '*************************************************************'
print

# Plotting out the chains
for ii in xrange(0, mc.ndim):
    print mc.pam_names[ii], chain_med[ii, 0], ' +\sigma ', chain_med[ii, 1], ' -\sigma ', chain_med[ii, 2]
    file_name = dir_output + 'chain_' + repr(ii) + '_' + mc.pam_names[ii] + '.png'
    plt.plot(chain_T[:, :, ii], '-', alpha=0.5)
    plt.axvline(mc.nburnin, c='r')
    plt.axhline(chain_med[ii, 0], c='k')
    plt.savefig(file_name, bbox_inches='tight')
    plt.close()

print
print '*************************************************************'
print

if args.mp != 'False':

    print 'MEGA plot'
    # plotting mega-corner plot

    megaplot_dat = np.zeros([np.size(flatchain,axis=0), np.size(flatchain, axis=1)+1])
    megaplot_med = np.zeros(np.size(flatchain, axis=1)+1)
    megaplot_dat[:, :-1] = flatchain[:,:]
    megaplot_dat[:, -1] = flatlnprob[:]
    megaplot_med[:-1] = chain_med[:,0]
    megaplot_med[-1] = lnprob_median
    labels = mc.pam_names
    labels.extend('Lnprob')
    fig = corner.corner(megaplot_dat[:, :], labels=labels, truths=megaplot_med)
    fig.savefig(dir_output + "ALL_corners.pdf", bbox_inches='tight')
    plt.close()

    print
    print '*************************************************************'
    print

if args.t != 'False':
    for nd in xrange(0, mc.ndim):  # (0,ndim):
        out_absc = np.arange(0, mc.nburnin, 1)
        out_lines = np.zeros(mc.nburnin)
        for ii in xrange(20, mc.nburnin):
            out_lines[ii] = GelmanRubin(chain_T[:ii, :, nd])

        #plt.ylim(0.95, 2.3)
        plt.plot(out_absc[20:], out_lines[20:], '-', color='k')
        plt.axhline(1.01)
        plt.savefig(dir_output + 'GRtrace_pam_' + repr(nd) + '.png', bbox_inches='tight')
        plt.close()

    print
    print '*************************************************************'
    print

if args.cc != 'False':
    cc = ChainConsumer()
    for nd in xrange(0, mc.ndim):  # (0,ndim):
        cc.add_chain(chain[:, :, nd].flatten(), walkers=mc.nwalkers)

    #print(cc.get_latex_table())
    print cc.get_summary()

    print cc.diagnostic_gelman_rubin(threshold=0.05)
    print cc.diagnostic_geweke()
    print
    print '*************************************************************'
    print

x0 = 1. / 150

M_star1_rand = np.random.normal(M_star1, M_star1_err, n_kept)

if 'kepler' in mc.model_list:

    if args.p != 'False' or args.v != 'False':

        plot_dir = dir_output + '/files_plot/'

        if not os.path.exists(plot_dir):
            os.makedirs(plot_dir)

        boundaries = np.asarray([mc.Tref, mc.Tref])
        plot_dir = dir_output + '/files_plot/'

        if not os.path.exists(plot_dir):
            os.makedirs(plot_dir)
        for dataset in mc.dataset_list:
            if dataset.kind == 'RV':
                boundaries[0] = min(boundaries[0], dataset.x[0])
                boundaries[1] = max(boundaries[1], dataset.x[-1])
        # tenpercent = (boundaries[1] - boundaries[0]) / 10.
        # boundaries = [boundaries[0] - tenpercent, boundaries[1] + tenpercent]
        boundaries += np.asarray([-1.0, 1.0])*(boundaries[1] - boundaries[0]) / 10.

        x_range = np.arange(boundaries[0], boundaries[1], 0.01)
        x_phase = np.arange(-0.50, 1.50, 0.005, dtype=np.double)

        model_dsys, model_plan, model_orbs, model_actv, model_curv = mc.rv_make_model(chain_med[:, 0], x_range, x_phase)

        fileout = open(plot_dir + 'curvature.dat', 'w')
        fileout.write('descriptor x_range m_curv \n')
        for ii in xrange(0, np.size(x_range)):
            fileout.write('{0:14f} {1:14f} \n'.format(x_range[ii], model_curv['BJD'][ii]))
        fileout.close()

    for planet_name in mc.pcv.planet_name:

        print 'Planet ', planet_name, ' summary'
        print

        dynamical_flag = (planet_name in mc.pcv.dynamical)


        """Let's feed the conversion function with the average results to get out a list of human variables"""
        convert_out = mc.pcv.convert(planet_name, chain_med)
        n_orbital = len(mc.pcv.var_list[planet_name])
        n_fitted = len(mc.pcv.var_list[planet_name]) - len(mc.pcv.fix_list[planet_name])

        n_curv = mc.ccv.order

        sample_plan = np.zeros([n_kept, n_orbital+6+n_curv])

        """Let's put all the human variables - including those that have been fixed - in sample_plan"""
        """An index is assigned to each variable to keep track of them in   """
        for n_var, var in enumerate(convert_out):
            convert_out[var] = n_var

        for ii in xrange(0, n_kept):

            convert_tmp = mc.pcv.convert(planet_name, flatchain[ii, :])
            for var in convert_out:
                sample_plan[ii, convert_out[var]] = convert_tmp[var]

        convert_out['Tperi'] = n_orbital + 1
        convert_out['Tcent'] = n_orbital + 2
        convert_out['M_kep'] = n_orbital + 3
        convert_out['a_smj'] = n_orbital + 4

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
            for n_var, var in enumerate(mc.ccv.list_pams):
                convert_out[var] = n_orbital + 6 + n_var

            for ii in xrange(0, n_kept):
                convert_tmp = mc.ccv.convert(flatchain[ii, :])
                for var in mc.ccv.list_pams:
                    sample_plan[ii, convert_out[var]] = convert_tmp[var]

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
            if mc.pcv.inclination[planet_name][1] > 0.01:
                sample_plan[:, convert_out['M_kep']] = sample_plan[:, convert_out['M_kep']] / np.sin(np.pi / 180. *
                        np.random.normal(mc.pcv.inclination[planet_name][0], mc.pcv.inclination[planet_name][1], n_kept))

        sample_med = np.asarray(map(lambda v: (v[1], v[2] - v[1], v[1] - v[0]),
                                zip(*np.percentile(sample_plan[:, :], [15.865, 50, 84.135], axis=0))))
        e_med = np.percentile(sample_plan[:, convert_out['e']], 65.865, axis=0)

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
            sel_label.append('P')
        if 'K' in mc.variable_list[planet_name]:
            sel_list.append(convert_out['K'])
            sel_label.append('K')

        if 'f' in mc.variable_list[planet_name]:
            sel_list.append(convert_out['f'])
            sel_label.append('phase')

        if 'e' in mc.variable_list[planet_name]:
            sel_list.append(convert_out['e'])
            sel_label.append('e')
        if 'o' in mc.variable_list[planet_name]:
            sel_list.append(convert_out['o'])
            sel_label.append('omega')
        if 'ecoso' in mc.variable_list[planet_name]:
            sel_list.append(convert_out['e'])
            sel_label.append('e')
            sel_list.append(convert_out['o'])
            sel_label.append('omega')

        if 'M' in mc.variable_list[planet_name]:
            sel_list.append(convert_out['M'])
            sel_label.append('M_e')
        else:
            sel_list.append(convert_out['M_kep'])
            sel_label.append('M_j')

        if 'curvature' in mc.model_list:
            for var in mc.ccv.list_pams:
                sel_list.append(convert_out[var])
                sel_label.append(var)

        fig = corner.corner(sample_plan[:, sel_list], labels=sel_label, truths=sample_med[sel_list, 0])
        fig.savefig(dir_output + planet_name + "_corners.pdf", bbox_inches='tight')
        plt.close()

        if args.p != 'False':
            # Write down the residuals
            color_list = ['b', 'g', 'r', 'c', 'm', 'y', 'k']

            f1, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=True)
            f2, (ax3, ax4) = plt.subplots(2, sharex=True, sharey=True)

            ax1.plot(x_range, model_plan['BJD'][planet_name], c='g')
            ax3.plot(x_phase, model_plan['pha'][planet_name], c='g')

            color_count = 0
            for dataset in mc.dataset_list:
                if dataset.kind == 'RV':
                    col_sel = color_list[color_count % 7]
                    color_count += 1
                    p_pha = (dataset.x0 / sample_med[0, 0]) % 1
                    y_det = dataset.y - model_dsys[dataset.name_ref] - model_curv[dataset.name_ref]
                    y_res = dataset.y - model_dsys[dataset.name_ref] - \
                            model_orbs[dataset.name_ref] - model_curv[dataset.name_ref]
                    y_1pl = y_res + model_plan[dataset.name_ref][planet_name]
                    ax1.errorbar(dataset.x, y_1pl, yerr=dataset.e, fmt=col_sel + '.', zorder=2)
                    #ax1.errorbar(dataset.x, y_det, yerr=dataset.e, fmt=col_sel + '.', zorder=2)
                    ax2.errorbar(dataset.x, y_res, yerr=dataset.e, fmt=col_sel + '.', zorder=2)

                    ax3.errorbar(p_pha, y_1pl, yerr=dataset.e, fmt=col_sel + '.', zorder=2)
                    ax4.errorbar(p_pha, y_res, yerr=dataset.e, fmt=col_sel + '.', zorder=2)

                    fileout = open(plot_dir + planet_name + '_' + dataset.name_ref + '_kep.dat', 'w')
                    fileout.write('descriptor BJD pha RV,+- RVdet,+- RVpla,+- RVres,+- RVmod,+- \n')
                    for ii in xrange(0, dataset.n):
                        fileout.write('{0:14f} {1:14f} {2:14f} {3:14f} {4:14f} {5:14f} {6:14f} {7:14f} '
                                      '{8:14f} {9:14f} {10:14f} {11:14f}'
                                      '\n'.format(dataset.x[ii], p_pha[ii],
                                                  dataset.y[ii], dataset.e[ii],
                                                  y_det[ii], dataset.e[ii], y_1pl[ii], dataset.e[ii],
                                                  y_res[ii], dataset.e[ii],
                                                  model_orbs[dataset.name_ref][ii], dataset.e[ii]))
                    fileout.close()

            f1.subplots_adjust(hspace=0)
            f1.savefig(plot_dir + planet_name + '_kep.pdf', bbox_inches='tight')
            f2.subplots_adjust(hspace=0)
            f2.savefig(plot_dir + planet_name + '_pha.pdf', bbox_inches='tight')
            plt.close()

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

        """ TO BE FIXED """
        """ Unsupported OLD section before using dictionaries """
        if 3 > 5:

            sample_plan = np.zeros([n_kept, 9])

            for ii in xrange(0, n_kept):
                convert_out = mc.pcv.convert(planet_name, flatchain[ii, :])
                sample_plan[ii, 0:5] = [convert_out['P'], convert_out['K'], convert_out['f'], convert_out['e'], convert_out['o']]

            # Time of periastron
            sample_plan[:, 6] = mc.Tref + (-sample_plan[:, 2] + sample_plan[:, 4]) / (2*np.pi) * sample_plan[:, 0]

            # Time of transit center
            sample_plan[:, 7] = mc.Tref + kp.kepler_Tcent_T0P(sample_plan[:, 0], sample_plan[:, 2], sample_plan[:, 3],
                                                                 sample_plan[:, 4])
            for ii in xrange(0, n_kept):
                # Planet mass
                sample_plan[ii, 5] = mc.M_SJratio * scipy.optimize.fsolve(get_mass, x0, args=(
                    M_star1_rand[ii], sample_plan[ii, 0], sample_plan[ii, 1], sample_plan[ii, 3]))
                # semi-major axis
                sample_plan[ii, 8] = np.power(
                    (Mu_sun * np.power(sample_plan[ii, 0] * seconds_in_day / (2 * np.pi), 2) / (AU_km ** 3.0)) *
                    M_star1_rand[ii], 1.00 / 3.00)

            if mc.pcv.inclination[planet_name][1] > 0.01:
                sample_plan[:, 5] = sample_plan[:, 5] / np.sin(np.pi/180. * np.random.normal(
                    mc.pcv.inclination[planet_name][0], mc.pcv.inclination[planet_name][1], n_kept))
                print ' Orbital inclination included in mass computation '

            sample_med = np.asarray(map(lambda v: (v[1], v[2] - v[1], v[1] - v[0]),
                                        zip(*np.percentile(sample_plan[:, :], [15.865, 50, 84.135], axis=0))))
            e_med = np.percentile(sample_plan[:, 3], 65.865, axis=0)

            '''
            sample_med[3, 0] = e_med

            sel_omega = (sample_plan[:, 4] > 0)
            omega_med = np.asarray(map(lambda v: (v[1], v[2] - v[1], v[1] - v[0]),
                                        zip(*np.percentile(sample_plan[sel_omega, :], [15.865, 50, 84.135], axis=0))))
            sample_med[4, :] = omega_med[4, :]

            esino_sel = flatchain[sel_omega, mc.variable_list[planet_name]['esino']]
            ecoso_sel = flatchain[sel_omega, mc.variable_list[planet_name]['ecoso']]

            print 'esino_sel (percentiles)', np.percentile(esino_sel, [15.865, 50, 84.135])
            print 'ecoso_sel (percentiles)', np.percentile(ecoso_sel, [15.865, 50, 84.135])
            '''

            print 'Period = ', sample_med[0, 0], ' +\sigma ', sample_med[0, 1], ' -\sigma ', sample_med[0, 2]
            print 'K      = ', sample_med[1, 0], ' +\sigma ', sample_med[1, 1], ' -\sigma ', sample_med[1, 2]
            print 'phase  = ', sample_med[2, 0], ' +\sigma ', sample_med[2, 1], ' -\sigma ', sample_med[2, 2]
            print 'e      = ', sample_med[3, 0], ' +\sigma ', sample_med[3, 1], ' -\sigma ', sample_med[3, 2], ', < ', e_med
            print 'o      = ', sample_med[4, 0], ' +\sigma ', sample_med[4, 1], ' -\sigma ', sample_med[4, 2]
            print 'Mass_J = ', sample_med[5, 0], ' +\sigma ', sample_med[5, 1], ' -\sigma ', sample_med[5, 2]
            print 'Mass_E = ', sample_med[5, 0]*mc.M_JEratio, \
                ' +\sigma ', sample_med[5, 1]*mc.M_JEratio, ' -\sigma ', sample_med[5, 2]*mc.M_JEratio
            print 'Tperi  = ', sample_med[6, 0], ' +\sigma ', sample_med[6, 1], ' -\sigma ', sample_med[6, 2]
            print 'Tcent  = ', sample_med[7, 0], ' +\sigma ', sample_med[7, 1], ' -\sigma ', sample_med[7, 2]
            print 'a      = ', sample_med[8, 0], ' +\sigma ', sample_med[8, 1], ' -\sigma ', sample_med[8, 2]


            sel_list = []
            sel_label = []

            #Left here for compatibility reasons
            if 'logP' in mc.variable_list[planet_name]:
                sel_list.append(0)
                sel_label.append('P')
            if 'logK' in mc.variable_list[planet_name]:
                sel_list.append(1)
                sel_label.append('K')
            if 'phase' in mc.variable_list[planet_name]:
                sel_list.append(2)
                sel_label.append('phase')

            if 'P' in mc.variable_list[planet_name]:
                sel_list.append(0)
                sel_label.append('P')
            if 'K' in mc.variable_list[planet_name]:
                sel_list.append(1)
                sel_label.append('K')
            if 'f' in mc.variable_list[planet_name]:
                sel_list.append(2)
                sel_label.append('phase')

            if 'e' in mc.variable_list[planet_name]:
                sel_list.append(3)
                sel_label.append('e')
            if 'o' in mc.variable_list[planet_name]:
                sel_list.append(4)
                sel_label.append('omega')
            if 'ecoso' in mc.variable_list[planet_name]:
                sel_list.append(3)
                sel_label.append('e')
                sel_list.append(4)
                sel_label.append('omega')

            sel_list.append(5)
            sel_label.append('Mass_J')

            #sel_list = [0, 1, 2, 3, 4, 6]
            #sel_label = ['P', 'K', 'phase', 'e', 'omega', 'Mass_J']

            fig = corner.corner(sample_plan[:, sel_list], labels=sel_label, truths=sample_med[sel_list, 0])
            fig.savefig(dir_output + planet_name + "_corners.pdf", bbox_inches='tight')
            plt.close()

            if args.p != 'False':
                # Write down the residuals
                color_list = ['b', 'g', 'r', 'c', 'm', 'y', 'k']

                f1, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=True)
                f2, (ax3, ax4) = plt.subplots(2, sharex=True, sharey=True)

                ax1.plot(x_range, model_plan['BJD'][planet_name], c='g')
                ax3.plot(x_phase, model_plan['pha'][planet_name], c='g')

                color_count = 0
                for dataset in mc.dataset_list:
                    if dataset.kind == 'RV':
                        col_sel = color_list[color_count % 7]
                        color_count += 1
                        p_pha = (dataset.x0 / sample_med[0, 0]) % 1
                        y_det = dataset.y-model_dsys[dataset.name_ref]
                        y_res = dataset.y-model_dsys[dataset.name_ref]-model_orbs[dataset.name_ref]
                        y_1pl = y_res + model_plan[dataset.name_ref][planet_name]

                        ax1.errorbar(dataset.x, y_1pl, yerr=dataset.e, fmt=col_sel+'.', zorder=2)
                        ax2.errorbar(dataset.x, y_res, yerr=dataset.e, fmt=col_sel+'.', zorder=2)

                        ax3.errorbar(p_pha, y_1pl, yerr=dataset.e, fmt=col_sel+'.', zorder=2)
                        ax4.errorbar(p_pha, y_res, yerr=dataset.e, fmt=col_sel+'.', zorder=2)

                        fileout = open(plot_dir + planet_name + '_' + dataset.name_ref + '_kep.dat','w')
                        fileout.write('descriptor BJD pha RV,+- RVdet,+- RVpla,+- RVres,+- RVmod,+- \n')
                        for ii in xrange(0, dataset.n):
                            fileout.write('{0:14f} {1:14f} {2:14f} {3:14f} {4:14f} {5:14f} {6:14f} {7:14f} '
                                          '{8:14f} {9:14f} {10:14f} {11:14f}'
                                          '\n'.format(dataset.x[ii], p_pha[ii],
                                                      dataset.y[ii], dataset.e[ii],
                                                      y_det[ii],  dataset.e[ii], y_1pl[ii], dataset.e[ii],
                                                      y_res[ii], dataset.e[ii],
                                                      model_orbs[dataset.name_ref][ii], dataset.e[ii]))
                        fileout.close()
                #'\n'.format(dataset.x[ii], (dataset.x0[ii] / sample_med[0, 0]) % 1,

                f1.subplots_adjust(hspace=0)
                f1.savefig(plot_dir + planet_name + '_kep.pdf', bbox_inches='tight')
                f2.subplots_adjust(hspace=0)
                f2.savefig(plot_dir + planet_name + '_pha.pdf', bbox_inches='tight')
                plt.close()
                #
                #f1, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=True)
                #ax1.plot(x_kep, model_plan['BJD'][planet_name], c='g')
                #for dataset in mc.dataset_list:
                #    if dataset.kind == 'RV':
                #        ax1.errorbar(dataset.x, rv_res+rv_pla, yerr=dataset.e, fmt='b.', zorder=2)
                #        ax2.errorbar(dataset.x, rv_res, yerr=dataset.e, fmt='b.', zorder=2)
                #f1.subplots_adjust(hspace=0)
                #f1.savefig(plot_dir + planet_name + '_kep.pdf', bbox_inches='tight')
                #plt.close()

                #f1, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=True)
                #ax1.plot(x_pha, y_pha, c='g')
                #for dataset in mc.dataset_list:
                #    if dataset.kind == 'RV':
                #        pp_pha = (dataset.x0 / sample_med[0, 0]) % 1
                #        ax1.errorbar(pp_pha, rv_res+rv_pla, yerr=dataset.e, fmt='b.', zorder=2)
                #        ax2.errorbar(pp_pha, rv_res, yerr=dataset.e, fmt='b.', zorder=2)
                #f1.subplots_adjust(hspace=0)
                #f1.savefig(plot_dir + planet_name + '_pha.pdf', bbox_inches='tight')
                #plt.close()

                fileout = open(plot_dir + planet_name + '_kep.dat', 'w')
                fileout.write('descriptor x_range m_keplr \n')
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
                if not os.path.exists(veusz_dir):
                    os.makedirs(veusz_dir)

                list_labels = ['P', 'K', 'e', 'm']

                n_int = 4

                output_plan = np.zeros([n_kept, n_int], dtype=np.double)
                output_plan[:, 0] = sample_plan[:, 0]
                output_plan[:, 1] = sample_plan[:, 1]
                output_plan[:, 2] = sample_plan[:, 4]
                output_plan[:, 3] = sample_plan[:, 5]
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
                            hist1d_norm = hist1d[0]*1. / n_kept
                            x_edges_1d = (x_edges[1:]+ x_edges[:-1])/2
                            data_grp.create_dataset(list_labels[ii]+'_x', data=x_edges_1d, compression="gzip")
                            data_grp.create_dataset(list_labels[ii]+'_y', data=hist1d_norm, compression="gzip")

                #plot lower-upper limits for the mass

                pams_limits = np.zeros([n_int, 2])
                for ii in xrange(0, n_int):
                    pams_limits[ii, :] = np.percentile(sample_plan[:, ii], [0.135, 99.865])

                random_n = 10000
                ii_kep = np.random.randint(low=0, high=n_kept, size=random_n)

                x_kep = np.arange(5900, 7900, 2)
                y_kep = np.zeros([np.size(x_kep), 3])

                x_pha = np.arange(-1.00, 2.00, 0.005,dtype=np.double)
                y_pha = np.zeros([np.size(x_pha), 3])

                y_flg = np.ones(random_n, dtype=bool)

                y_kep_tmp = np.zeros(np.size(x_kep))
                y_pha_tmp = np.zeros(np.size(x_pha))

                if pams_limits[3, 0] < 0.02 : pams_limits[3, 0] = 0.00
                for ii in xrange(0, n_int):
                    y_flg = y_flg & (sample_plan[ii_kep, ii] >= pams_limits[ii,0]) & (sample_plan[ii_kep, ii] <= pams_limits[ii, 1])

                y_kep[:, 2] = kp.kepler_RV_T0P(x_kep - mc.Tref, sample_med[2, 0], sample_med[0, 0], sample_med[1, 0],
                                               sample_med[3, 0], sample_med[4, 0])
                y_kep[:, 0] = y_kep[:, 2]
                y_kep[:, 1] = y_kep[:, 2]
                # value initialization
                y_pha[:, 2] = kp.kepler_RV_T0P(x_pha*sample_med[0, 0], sample_med[2, 0], sample_med[0, 0], sample_med[1, 0],
                                               sample_med[3, 0], sample_med[4, 0])
                y_pha[:, 0] = y_pha[:, 2]
                y_pha[:, 1] = y_pha[:, 2]

                print 'Accepted randomizations: ', np.sum(y_flg)
                for ii in xrange(0, random_n):
                    ir = ii_kep[ii]
                    if y_flg[ii]:
                        y_kep_tmp = kp.kepler_RV_T0P(x_kep - mc.Tref, sample_plan[ir, 2], sample_plan[ir, 0], sample_plan[ir, 1], sample_plan[ir, 3], sample_plan[ir, 4])
                        y_pha_tmp = kp.kepler_RV_T0P(x_pha*sample_plan[ir, 0], sample_plan[ir, 2], sample_plan[ir, 0], sample_plan[ir, 1], sample_plan[ir, 3], sample_plan[ir, 4])

                        y_kep[:, 0] = np.minimum(y_kep[:, 0],y_kep_tmp)
                        y_kep[:, 1] = np.maximum(y_kep[:, 1],y_kep_tmp)

                        y_pha[:, 0] = np.minimum(y_pha[:, 0],y_pha_tmp)
                        y_pha[:, 1] = np.maximum(y_pha[:, 1],y_pha_tmp)

                fig1 = plt.plot(x_kep, y_kep[:, 2], c='g')
                fig1 = plt.plot(x_kep, y_kep[:, 0], c='r')
                fig1 = plt.plot(x_kep, y_kep[:, 1], c='b')
                plt.savefig(veusz_dir + planet_name + '_kep.png', bbox_inches='tight')
                plt.close()

                fig2 = plt.plot(x_pha, y_pha[:, 2], c='g')
                fig2 = plt.plot(x_pha, y_pha[:, 0], c='r')
                fig2 = plt.plot(x_pha, y_pha[:, 1], c='b')
                plt.savefig(veusz_dir + planet_name + '_pha.png', bbox_inches='tight')
                plt.close()

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

if 'gaussian' in mc.model_list:
    n_vars = 0
    sample_plan_transpose = []
    sel_label = []

    for name in mc.gcv.list_pams_common:
        if name in mc.variable_list['Common']:
            n_vars += 1
            # mc.gcv.var_list[name] select the corresponding value in the emcee input array (theta)
            var = flatchain[:, mc.gcv.var_list['Common'][name]]
            # mc.gcv.variables[name](theta, mc.gcv.fixed, :) convert the value into physical unit
            # e.g. from logarithmic to linear space)
            var_phys = mc.gcv.variables['Common'][name](var, var, xrange(0, n_kept))
            sample_plan_transpose.append(var_phys)
            sel_label.append(name)

    for dataset in mc.dataset_list:
        for name in mc.gcv.list_pams_dataset:
            n_vars += 1
            var = flatchain[:, mc.gcv.var_list[dataset.name_ref][name]]
            var_phys = mc.gcv.variables[dataset.name_ref][name](var, var, xrange(0, n_kept))
            sample_plan_transpose.append(var_phys)
            sel_label.append(dataset.name_ref + '_' + name)

    sample_plan = np.asarray(sample_plan_transpose).T
    sample_med = np.asarray(map(lambda v: (v[1], v[2] - v[1], v[1] - v[0]),
                                    zip(*np.percentile(sample_plan[:, :], [15.865, 50, 84.135], axis=0))))

    for lab, sam in zip(sel_label, sample_med):
        print lab,' = ', sam[0], ' +\sigma ', sam[1], ' -\sigma ', sam[2]
    print

    fig = corner.corner(sample_plan[:, :], labels=sel_label, truths=sample_med[:, 0])
    fig.savefig(dir_output + "GPs_corners.pdf", bbox_inches='tight')
    plt.close()


if 'curvature' in mc.model_list:
    n_vars = 0
    sample_plan_transpose = []
    sel_label = []

    for name in mc.ccv.list_pams:
        n_vars += 1
        var = flatchain[:, mc.ccv.var_list[name]]
        var_phys = mc.ccv.variables[name](var, var, xrange(0, n_kept))
        sample_plan_transpose.append(var_phys)
        sel_label.append(name)

    sample_plan = np.asarray(sample_plan_transpose).T
    sample_med = np.asarray(map(lambda v: (v[1], v[2] - v[1], v[1] - v[0]),
                                    zip(*np.percentile(sample_plan[:, :], [15.865, 50, 84.135], axis=0))))

    for lab, sam in zip(sel_label, sample_med):
        print lab,' = ', sam[0], ' +\sigma ', sam[1], ' -\sigma ', sam[2]
    print

print
print 'chain shape ', np.shape(chain)