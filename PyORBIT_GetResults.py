from classes.model_container import ModelContainer
from classes.input_parser import yaml_parser, pars_input
from classes.io_subroutines import pyde_save_to_pickle, pyde_load_from_cpickle, \
    emcee_save_to_cpickle, emcee_load_from_cpickle, emcee_flatchain, emcee_flatlnprob, \
    GelmanRubin, model_container_plot
import numpy as np
import emcee
from pyde.de import DiffEvol
import os
import argparse
import matplotlib as mpl
from matplotlib.ticker import FormatStrFormatter
import sys
mpl.use('Agg')
from matplotlib import pyplot as plt
import corner
import classes.constants as constants


def pyorbit_getresults(config_in, sampler, plot_dictionary):

    #plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman']})
    plt.rcParams["font.family"] = "Times New Roman"
    plt.rc('text', usetex=True)

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


    sample_keyword = {
        'polychord':['polychord', 'PolyChord', 'polychrod', 'poly'],
        'emcee': ['emcee', 'MCMC', 'Emcee']
    }

    if sampler in sample_keyword['emcee']:

        dir_input = './' + config_in['output'] + '/emcee/'
        dir_output = './' + config_in['output'] + '/emcee_plot/'
        os.system('mkdir -p ' + dir_output)

        mc, starting_point, population, prob, state, \
        sampler_chain, sampler_lnprobability, sampler_acceptance_fraction = \
            emcee_load_from_cpickle(dir_input)

        mc.model_setup()
        """ Required to create the right objects inside each class - if defined inside """
        theta_dictionary = mc.get_theta_dictionary()

        nburnin = mc.emcee_parameters['nburn']
        nthin = mc.emcee_parameters['thin']

        flat_chain = emcee_flatchain(sampler_chain, nburnin, nthin)
        flat_lnprob = emcee_flatlnprob(sampler_lnprobability, nburnin, nthin)

        chain_med = np.asarray(map(lambda v: (v[1], v[2] - v[1], v[1] - v[0]),
                                   zip(*np.percentile(flat_chain[:, :], [15.865, 50, 84.135], axis=0))))
        lnprob_med = np.percentile(flat_lnprob, [15.865, 50, 84.135], axis=0)
        lnprob_med[1:] = np.abs(lnprob_med[1:] - lnprob_med[0])

        print
        print 'Reference Time Tref: ', mc.Tref
        print
        print 'Dimensions = ', mc.ndim
        print 'Nwalkers = ', mc.emcee_parameters['nwalkers']
        print
        print '*************************************************************'
        print

        lnprob_median = np.median(flat_lnprob)

        print 'median log-likelihood: ', lnprob_median
        #mc.results_resumen(flat_chain)

        fig = plt.figure(figsize=(12, 12))
        plt.xlabel('$\ln \mathcal{L}$')
        plt.plot(sampler_lnprobability.T, '-', alpha=0.5)
        plt.axhline(lnprob_median)
        plt.axvline(nburnin/nthin, c='r')
        plt.savefig(dir_output + 'LNprob_chain.png', bbox_inches='tight', dpi=300)
        plt.close(fig)

        print
        print '*************************************************************'
        print

        if plot_dictionary['full_correlation']:

            print 'full_correlation plot'
            # plotting mega-corner plot
            plt.rc('text', usetex=False)

            full_corner_dat = np.zeros([np.size(flat_chain, axis=0), np.size(flat_chain, axis=1) + 1])
            full_corner_med = np.zeros(np.size(flat_chain, axis=1) + 1)
            full_corner_dat[:, :-1] = flat_chain[:, :]
            full_corner_dat[:, -1] = flat_lnprob[:]
            full_corner_med[:-1] = chain_med[:, 0]
            full_corner_med[-1] = lnprob_median

            full_corner_labels = [repr(ii) for ii in xrange(0, np.size(flat_chain, axis=0))]

            for theta_name, ii in theta_dictionary.iteritems():
                full_corner_labels[ii] = theta_name
            full_corner_labels.extend(['ln L'])
            fig = corner.corner(full_corner_dat[:, :], labels=full_corner_labels, truths=full_corner_med)
            fig.savefig(dir_output + "full_corners.pdf", bbox_inches='tight', dpi=300)
            plt.close(fig)
            plt.rc('text', usetex=True)

            print
            print '*************************************************************'
            print

        if plot_dictionary['chains']:
            print 'plotting the chains... '
            for theta_name, ii in theta_dictionary.iteritems():
                file_name = dir_output + 'chain_' + repr(ii) + '_' + theta_name + '.png'
                fig = plt.figure(figsize=(12, 12))
                plt.plot(sampler_chain[:, :, ii].T, '-', alpha=0.5)
                plt.axvline(nburnin/nthin, c='r')
                plt.savefig(file_name, bbox_inches='tight', dpi=300)
                plt.close(fig)

            print
            print '*************************************************************'
            print

        if plot_dictionary['traces']:
            print 'plotting the traces... '
            for theta_name, th in theta_dictionary.iteritems():
                file_name = dir_output + 'GRtrace_' + repr(th) + '_' + theta_name + '.png'
                out_absc = np.arange(0, nburnin/nthin, 1)
                out_lines = np.zeros(nburnin/nthin)
                for ii in xrange(20, nburnin/nthin):
                    out_lines[ii] = GelmanRubin(sampler_chain[:, :ii, th].T)
                fig = plt.figure(figsize=(12, 12))
                plt.plot(out_absc[20:], out_lines[20:], '-', color='k')
                plt.axhline(1.01)
                plt.savefig(file_name, bbox_inches='tight', dpi=300)
                plt.close(fig)

            print
            print '*************************************************************'
            print


        print

        model_out, logchi2_out = mc.get_model(chain_med[:, 0])

        mc_deepcopy = model_container_plot(mc)
        model_plot, _ = mc_deepcopy.get_model(chain_med[:, 0])

        fig = plt.figure(figsize=(12, 12))
        for dataset_name, dataset in mc.dataset_dict.iteritems():
            plt.scatter(dataset.x0, dataset.y - model_out[dataset_name]['systematics'])

        for dataset_name, dataset in mc_deepcopy.dataset_dict.iteritems():
            plt.plot(dataset.x0, model_plot[dataset_name]['complete'], c='r')

        plt.savefig(dir_output + '_test.png', bbox_inches='tight', dpi=300)
        plt.close(fig)
        print

if __name__ == '__main__':
    print 'This program is being run by itself'

    parser = argparse.ArgumentParser(prog='PyORBIT_GetResults.py', description='PyDE+emcee runner')
    parser.add_argument('sample', type=str, nargs=1, help='sample (emcee or polychord)')
    parser.add_argument('config_file', type=str, nargs=1, help='config file')
    parser.add_argument('-p', type=str, nargs='?', default=False, help='Create plot files')
    parser.add_argument('-c', type=str, nargs='?', default=False, help='Create chains plots')
    parser.add_argument('-t', type=str, nargs='?', default=False, help='Create Gelman-Rubin traces')
    parser.add_argument('-m', type=str, nargs='?', default=False, help='Create full corellation plot - it may be slow!')

    plot_dictionary = {
        'chains': False,
        'plot': False,
        'traces': False,
        'full_correlation': False
    }

    args = parser.parse_args()
    sampler = args.sample[0]
    file_conf = args.config_file[0]

    if args.p is not False :
        plot_dictionary['plot'] = True
    if args.c is not False :
        plot_dictionary['chains'] = True
    if args.t is not False:
        plot_dictionary['traces'] = True
    if args.m is not False:
        plot_dictionary['full_correlation'] = True

    print plot_dictionary
    config_in = yaml_parser(file_conf)

    pyorbit_getresults(config_in, sampler, plot_dictionary)

else:
    print 'I am being imported from another module'
